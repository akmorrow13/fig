/**
 * Copyright 2015 Frank Austin Nothaft
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package net.fnothaft.fig.models

import net.fnothaft.fig.avro._
import org.apache.spark.SparkContext._
import org.apache.spark.Logging
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.ReferenceRegion._
import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import org.bdgenomics.adam.rdd.BroadcastRegionJoin
import org.bdgenomics.adam.rich.RichGenotype._
import org.bdgenomics.formats.avro.{ Genotype, GenotypeAllele, Variant }
import scala.annotation.tailrec
import scala.collection.JavaConversions._
import scala.collection.mutable.Buffer

private[models] case class Shift(location: Long,
                                 shift: Int) {
}

object VariantPromoter extends Serializable with Logging {

  def apply(pRdd: RDD[ReferencePromoter],
            gRdd: RDD[Genotype],
            motifRepository: MotifRepository): RDD[VariantPromoter] = {
    // cache promoter RDD; we'll need it a few times
    pRdd.cache()

    // region join genotypes versus promoters
    val jRdd = BroadcastRegionJoin.partitionAndJoin(pRdd.keyBy(_.pos),
                                                    gRdd.keyBy(ReferencePosition(_)
                                                      .asInstanceOf[ReferenceRegion]))

    // map down and group-by sample and gene
    // then join against original reference promoter RDD
    val sampleGeneGenotypes = jRdd.map(kv => {
      val (p, g) = kv

      ((p.geneId, g.getSampleId), g)
    }).groupByKey()
      .map(kv => {
        val ((gene, sample), genotypes) = kv
        (gene, (sample, genotypes))
      })
      .join(pRdd.keyBy(_.geneId))
   
    // unpersist promoter RDD
    pRdd.unpersist()

    // map to create variant promoters
    sampleGeneGenotypes.flatMap(kv => {
      val (geneId, ((sample, genotypes), refPromoter)) = kv
      VariantPromoter(geneId, sample, genotypes, refPromoter, motifRepository)
    })
  }

  private[models] def apply(geneId: String,
                            sampleId: String,
                            genotypes: Iterable[Genotype],
                            refPromoter: ReferencePromoter,
                            motifRepository: MotifRepository): Iterable[VariantPromoter] = {
    assert(genotypes.forall(_.getSampleId == sampleId),
           "Somehow, samples got mixed up in %s/%s.".format(geneId, sampleId))

    // get ploidy
    val ploidies = genotypes.map(g => g.getAlleles.size)
      .toSeq
      .distinct
    
    // get phasing info
    val isPhased = genotypes.forall(g => g.getIsPhased)

    if (ploidies.length != 1) {
      log.warn("Discarding genotypes on %s for sample %s, as variable ploidy was seen.".format(
        geneId,
        sampleId))
      Iterable.empty
    } else if (!isPhased && ploidies.head != 1) {
      log.warn("Discarding genotypes on %s for sample %s, as genotypes are unphased.".format(
        geneId,
        sampleId))
      Iterable.empty
    } else {
      // what is our ploidy?
      val ploidy = ploidies.head

      // per copy, filter out reference calls
      val variantsPerHaplotype = (0 until ploidy).map(i => {
        (i, genotypes.flatMap(g => {
          // we only keep this variant if it is the true alt call
          if (g.getAlleles.get(i) == GenotypeAllele.Alt) {
            Some(g.getVariant)
          } else {
            None
          }
        }))
      })

      try {
        variantsPerHaplotype.map(kv => {
          val (copy, variants) = kv
          VariantPromoter(refPromoter,
                          sampleId,
                          copy,
                          ploidy,
                          variants,
                          motifRepository)
        })
      } catch {
        case iae: IllegalArgumentException => {
          Iterable.empty
        }
      }
    }
  }

  private[models] def apply(promoter: ReferencePromoter,
                            sampleId: String,
                            copy: Int,
                            samplePloidy: Int,
                            variants: Iterable[Variant],
                            motifRepository: MotifRepository): VariantPromoter = {
    if (variants.isEmpty) {
      VariantPromoter(promoter,
                      sampleId,
                      copy,
                      samplePloidy,
                      promoter.sequence,
                      variants,
                      promoter.tfbs)
    } else {
      // flatten variants into haplotype
      val (sequence, shifts) = extractHaplotype(promoter.sequence,
                                                promoter.pos,
                                                variants)
      
      // reannotate tfbs
      val tfbs = annotateSites(sequence,
                               shifts,
                               promoter.tfbs,
                               promoter.pos,
                               motifRepository)

      VariantPromoter(promoter,
                      sampleId,
                      copy,
                      samplePloidy,
                      sequence,
                      variants,
                      tfbs)      
    }
  }

  private[models] def annotateSites(sequence: String,
                                    shifts: Iterable[Shift],
                                    sites: Iterable[BindingSite],
                                    region: ReferenceRegion,
                                    motifRepository: MotifRepository): Iterable[BindingSite] = {

    def getShift(loc: Long): Int = {
      shifts.filter(_.location <= loc)
        .map(_.shift)
        .fold(0)(_ + _)
    }
    
    sites.flatMap(site => {
      // does this site need to be shifted? if so, how far
      val shiftBy = getShift(site.getStart)
      
      // get the substring for this site
      val siteString = sequence.substring((site.getStart - region.start).toInt + shiftBy,
                                          (site.getEnd - region.start).toInt + shiftBy)

      // annotate site
      motifRepository.annotate(siteString,
                               site)
                                 .map(site => {
                                   site.setShift(shiftBy)
                                   site
                                 })
    })
  }

  private[models] def extractHaplotype(refSequence: String,
                                       refRegion: ReferenceRegion,
                                       variants: Iterable[Variant]): (String,
                                                                      Iterable[Shift]) = {

    // create a list to append shifts to
    var shiftList = List.empty[Shift]

    // string to iter
    val refIter = refSequence.toIterator

    // initialize our position with the head of our region
    var pos = refRegion.start

    // create a string builder
    val sb = new StringBuilder()

    // helper function for advancing by popping characters off of the iterator and
    // adding them to the string builder
    @tailrec def advance(iteration: Long) {
      if (iteration > 0L) {
        sb.append(refIter.next)
        advance(iteration - 1L)
      }
    }
    
    // helper method for moving to the new variant loci
    def moveTo(variant: Variant) {
      assert(variant.getContig.getContigName == refRegion.referenceName)
      assert(variant.getStart <= refRegion.end && variant.getStart > pos)

      // how far do we need to advance?
      val distance = variant.getStart - pos
      pos = variant.getEnd
      
      // advance
      advance(distance)
    }

    // sort variants, and then convert to an iterator
    val vIter = variants.toSeq
      .sortBy(ReferencePosition(_).asInstanceOf[ReferenceRegion])
      .toIterator

    // helper function for reconstrucing the ref sequence at a variant locus
    @tailrec def take(v: Variant,
                      i: Long = 0L,
                      vsb: StringBuilder = new StringBuilder()): String = {
      if (i >= (v.getEnd - v.getStart)) {
        vsb.toString
      } else {
        vsb.append(refIter.next)
        take(v, i + 1, vsb)
      }
    }

    @tailrec def crawl {
      // do we still have variants?
      if (vIter.hasNext) {
        val v = vIter.next

        // reposition to this variant
        moveTo(v)

        // the variant's ref must match the ref here, so let's check that
        val refSeq = take(v)
        assert(refSeq == v.getReferenceAllele,
               "Computed reference sequence (%s) did not match variant reference sequence (%s)."
                 .format(refSeq, v.getReferenceAllele))

        // append the alt allele
        sb.append(v.getAlternateAllele)

        // if the ref and alt are different lengths, we must add a shift
        val shift = v.getAlternateAllele.length - v.getReferenceAllele.length
        if (shift != 0) {
          shiftList = Shift(pos, shift) :: shiftList
        }
        
        // keep crawling
        crawl
      }
    }

    // start crawling
    crawl

    // helper function for advancing to the end of the iterator by popping characters
    // off of the iterator and adding them to the string builder
    @tailrec def advanceToEnd {
      if (refIter.hasNext) {
        sb.append(refIter.next)
        advanceToEnd
      }
    }

    // when we're done, pad out the string builder until we hit the end of the iterator
    advanceToEnd
    
    (sb.toString, shiftList.toIterable)
  }
}

case class VariantPromoter(promoter: ReferencePromoter,
                           sampleId: String,
                           copy: Int,
                           samplePloidy: Int,
                           sequence: String,
                           variants: Iterable[Variant],
                           tfbs: Iterable[BindingSite]) {

  val tfbsSpacing = {
    val tfbsSeq = tfbs.toSeq
    (0 until tfbs.size).flatMap(i => {
      ((i + 1) until tfbs.size).flatMap(j => {
        if (tfbsSeq(i).getEnd + tfbsSeq(i).getShift <=
          tfbsSeq(j).getStart + tfbsSeq(j).getShift) {
          Some((tfbsSeq(j).getStart + tfbsSeq(j).getShift) -
               (tfbsSeq(i).getEnd + tfbsSeq(i).getShift))
        } else if (tfbsSeq(i).getStart + tfbsSeq(i).getShift >=
          tfbsSeq(j).getEnd + tfbsSeq(j).getShift) {
          Some((tfbsSeq(i).getStart + tfbsSeq(i).getShift) -
               (tfbsSeq(j).getEnd + tfbsSeq(j).getShift))
        } else {
          None
        }
      })
    })
  }

  def label: LabeledPromoter = {
    // create builder for promoter
    val lpb = LabeledPromoter.newBuilder()
      .setSampleId(sampleId)
      .setGene(promoter.geneId)
      .setHaplotypeNumber(copy)
      .setCopyNumber(samplePloidy)
      .setVariants(bufferAsJavaList(variants.toBuffer))
      .setGcRatio(sequence.count(c => c == 'C' || c == 'G').toDouble /
                  sequence.length.toDouble)
      .setOriginalGcRatio(promoter.gcRatio)
      .setSpacingChanges(tfbsSpacing.count(s => s % 5 == 0) -
                         promoter.tfbsSpacing.count(s => s % 5 == 0))

    // now, let's loop over all tfbs and see which we lost or modified
    val keptSites = Buffer.empty[BindingSite]
    val modifiedSites = Buffer.empty[ModifiedBindingSite]
    val varSitesSeq = tfbs.toSeq
    var refSitesSeq = promoter.tfbs.toSeq

    varSitesSeq.foreach(vs => {
      def innerLoop {
        (0 until refSitesSeq.length).foreach(i => {
          val rs = refSitesSeq(i)
          
          // are we looking at the same site?
          if (vs.getStart == rs.getStart &&
              vs.getEnd == rs.getEnd &&
              vs.getOrientation == rs.getOrientation &&
              vs.getTf == rs.getTf) {
            
            // was this site modified?
            if (vs.getSequence == rs.getSequence) {
              keptSites += vs
            } else {
              modifiedSites += ModifiedBindingSite.newBuilder()
                .setTf(vs.getTf)
                .setContig(vs.getContig)
                .setStart(vs.getStart)
                .setEnd(vs.getEnd)
                .setShift(vs.getShift)
                .setOrientation(vs.getOrientation)
                .setAffinityChange(vs.getPredictedAffinity / rs.getPredictedAffinity)
                .build()
            }
            
            // filter this element out of the seq
            refSitesSeq = refSitesSeq.take(i) ++ refSitesSeq.drop(i + 1)
            
            // break the loop
            return
          }

        })
      }

      innerLoop
    })

    // add binding sites
    lpb.setUnmodifiedTfbs(bufferAsJavaList(keptSites))
      .setModifiedTfbs(bufferAsJavaList(modifiedSites))
      .setLostTfbs(bufferAsJavaList(refSitesSeq.toBuffer))

    // build and return
    lpb.build()
  }
}
