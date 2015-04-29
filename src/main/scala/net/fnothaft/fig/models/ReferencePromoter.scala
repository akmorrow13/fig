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

import net.fnothaft.fig.avro.BindingSite
import org.apache.spark.SparkContext._
import org.apache.spark.rdd.RDD
import org.bdgenomics.adam.models.{ Gene, ReferencePosition, ReferenceRegion }
import org.bdgenomics.adam.rdd.BroadcastRegionJoin
import org.bdgenomics.adam.util.ReferenceFile
import org.bdgenomics.formats.avro.{ Contig, Strand }
import scala.math.max

object ReferencePromoter extends Serializable {

  private[models] def extractBindingSites(rdd: RDD[String]): RDD[(ReferenceRegion, BindingSite)] = {
    // per feature, we want to:
    // - key each feature with the TF name
    // - extract start/end position
    rdd.map(l => {
      val sl = l.split(' ')
      val ctg = sl(1)
      val start = sl(2).toLong
      // end is inclusive
      val end = sl(3).toLong + 1L
      val orientation = if(sl(4) == "-") {
        Strand.Reverse
      } else if (sl(4) == "+") {
        Strand.Forward
      } else {
        throw new IllegalArgumentException("Unknown strand '%s' for: %s.".format(sl(4), l))
      }
      
      (ReferenceRegion(ctg,
                       start,
                       end),
       BindingSite.newBuilder()
         .setTf(sl(0).takeWhile(_ != '_'))
         .setContig(Contig.newBuilder().setContigName(ctg).build())
         .setStart(start)
         .setEnd(end)
         .setOrientation(orientation)
         .build())
    })
  }

  def apply(gRdd: RDD[Gene],
            features: RDD[String],
            twoBit: ReferenceFile,
            startDistance: Int,
            stopDistance: Int,
            motifRepository: MotifRepository): RDD[ReferencePromoter] = {
    require(startDistance > stopDistance,
            "The start distance (%d) must be greater than the stop distance (%d).".format(
      startDistance,
      stopDistance))

    // from gene rdd, we want to extract:
    // - the gene name
    // - the regulatory region before the gene
    val namePromoterTss = gRdd.map(g => {
      // for simplicity's sake, we're going to take the transcription start site
      // from the first transcript in the gene
      val tRegion = g.transcripts.head.region
      val tss = ReferencePosition(tRegion.referenceName, tRegion.start)

      // we add flanks to the TSS to get the promoter region
      val pRegion = ReferenceRegion(tss.referenceName,
                                    max(tss.start - startDistance, 0L),
                                    max(tss.start - stopDistance, 0L))

      // this will be used in a region join later, so key with the region
      (pRegion, (g.id, tss))
    })

    // we will cache the last RDD, as we'll use it in a region join _and_ a join
    namePromoterTss.cache()

    val bindingSites = extractBindingSites(features)

    // now, we run a region join. The genes have reasonably low cardinality (~30k),
    // so we will use a broadcast region join.
    val regionJoinedRdd = BroadcastRegionJoin.partitionAndJoin(namePromoterTss,
                                                               bindingSites)

    // we now need to group the binding sites by gene
    val sitesByGene = regionJoinedRdd.map(kvp => {
      val ((geneName, _), tfbs) = kvp
      (geneName, tfbs)
    }).groupByKey()

    // then, join the sites against the promoter RDD
    val joinedRdd = namePromoterTss.map(kvp => (kvp._2._1, (kvp._1, kvp._2._2)))
      .join(sitesByGene)

    // we can now unpersist the namePromoterTss RDD
    namePromoterTss.unpersist()

    // broadcast the two bit file
    val broadcast2Bit = joinedRdd.context.broadcast(twoBit)

    // finish by mapping to create the reference promoter RDD
    joinedRdd.map(kvp => {
      val (gene, ((region, position), tfbs)) = kvp

      val refSeq = broadcast2Bit.value.extract(region)
      val annotatedTfbs = motifRepository.annotate(tfbs, refSeq, region)

      ReferencePromoter(gene,
                        region,
                        position,
                        refSeq,
                        annotatedTfbs)
    })
  }
}

case class ReferencePromoter(geneId: String,
                             pos: ReferenceRegion,
                             tss: ReferencePosition,
                             sequence: String,
                             tfbs: Iterable[BindingSite]) {

  val gcRatio = (sequence.count(c => c == 'C' || c == 'G' || c == 'c' || c == 'g').toDouble /
                 sequence.length.toDouble)
  val tfbsSpacing = {
    val tfbsSeq = tfbs.toSeq
    (0 until tfbs.size).flatMap(i => {
      ((i + 1) until tfbs.size).flatMap(j => {
        if (tfbsSeq(i).getEnd <= tfbsSeq(j).getStart) {
          Some(tfbsSeq(j).getStart - tfbsSeq(i).getEnd)
        } else if (tfbsSeq(i).getStart >= tfbsSeq(j).getEnd) {
          Some(tfbsSeq(i).getStart - tfbsSeq(j).getEnd)
        } else {
          None
        }
      })
    })
  }
}
