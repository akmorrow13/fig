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
import org.apache.spark.{ Logging, SparkContext }
import org.apache.spark.broadcast.Broadcast
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.formats.avro.Strand

object MotifRepository extends Serializable {

  def apply(sc: SparkContext,
            motifs: Seq[Motif]): MotifRepository = {
    BroadcastMotifRepository(sc.broadcast(motifs.groupBy(_.label)))
  }
}

trait MotifRepository extends Serializable with Logging {

  val motifs: Map[String, Seq[Motif]]

  def annotate(tfbs: Iterable[BindingSite],
               sequence: String,
               region: ReferenceRegion): Iterable[BindingSite] = {
    tfbs.flatMap(site => {
      assert(site.getContigName == region.referenceName)

      // is this site fully contained in the region?
      if (site.getStart >= region.start &&
          site.getEnd <= region.end) {
        val motifSequence = sequence.substring((site.getStart - region.start).toInt,
                                               (site.getEnd - region.start).toInt)
        annotate(motifSequence, site)
      } else {
        // print log message and drop
        log.warn("Dropping TFBS (%s) that falls outside of region %s.".format(site,
                                                                              region))
        None
      }
    })
  }           

  def annotate(sequence: String,
               tfbs: BindingSite): Option[BindingSite] = {
    require(tfbs.getEnd - tfbs.getStart == sequence.length,
            "Sequence length (%s) is incompatible with binding site (%s).".format(sequence, tfbs))

    // find motif class
    val motifClass = motifs.getOrElse(tfbs.getTf,
                                      throw new IllegalArgumentException("Could not find motif for %s.".format(tfbs.getTf)))

    // search for correct length motif
    val correctLengthMotifs = motifClass.filter(_.length == sequence.length)
    require(correctLengthMotifs.length != 0,
            "Did not find a correct length motif for %s with sequence length %d.".format(tfbs.getTf,
                                                                                         sequence.length))
    require(correctLengthMotifs.length == 1,
            "Found multiple motifs (%s) with sequence length %d for %s.".format(correctLengthMotifs.mkString(","),
                                                                                sequence.length,
                                                                                tfbs.getTf))
    val motif = correctLengthMotifs.head

    // score the motif
    val (motifScore, motifSequence) = tfbs.getOrientation match {
      case Strand.FORWARD => (motif.sequenceProbability(sequence), sequence)
      case Strand.REVERSE => {
        val rev = sequence.reverse
        (motif.sequenceProbability(rev), rev)
      }
    }

    // did we get a non-zero score? if so, annotate the tfbs, else return an empty option
    if (motifScore != 0.0) {
      Some(BindingSite.newBuilder(tfbs)
        .setSequence(motifSequence)
        .setPredictedAffinity(motifScore)
        .build())
    } else {
      None
    }
  }
}

case class BroadcastMotifRepository(private val bcastMotifsByName: Broadcast[Map[String, Seq[Motif]]]) extends MotifRepository {

  lazy val motifs = bcastMotifsByName.value
}
