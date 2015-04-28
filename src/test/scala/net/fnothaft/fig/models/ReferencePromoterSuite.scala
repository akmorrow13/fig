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

import net.fnothaft.fig.FigFunSuite
import org.apache.spark.SparkContext
import org.bdgenomics.adam.models._
import org.bdgenomics.adam.util.ReferenceFile
import org.bdgenomics.formats.avro.{ Contig, Feature, Strand }
import org.bdgenomics.utils.misc.MathUtils

case class SimpleReferenceFile() extends ReferenceFile {
  //         ABCDEFGHIJKLMNOPQRSTUVWXYZ
  val ref = "AACTGCGACCTATCTGGGATACTCGC"
  def extract(region: ReferenceRegion): String = {
    require(region.referenceName == "alphabet")
    ref.substring(region.start.toInt, region.end.toInt)
  }
}

class ReferencePromoterSuite extends FigFunSuite {

  val rf = SimpleReferenceFile()

  def makeAnnotator(sc: SparkContext): MotifRepository = {
    val motifs = Seq(Motif("tf1", Array(0.0, 1.0, 0.0, 0.0, // C
                                        0.0, 0.0, 0.0, 1.0, // T
                                        0.0, 0.0, 1.0, 0.0) // G
                           ),
                     Motif("tf2", Array(0.0, 0.0, 0.0, 1.0, // T
                                        0.5, 0.5, 0.0, 0.0, // M
                                        0.0, 0.0, 0.0, 1.0) // T
                           ),
                     Motif("tf3", Array(0.5, 0.5, 0.0, 0.0, // M
                                        0.5, 0.5, 0.0, 0.0, // M
                                        0.5, 0.5, 0.0, 0.0) // M
                           ))
    MotifRepository(sc, motifs)
  }

  test("test extraction from test file") {
    intercept[IllegalArgumentException] {
      rf.extract(ReferenceRegion("1", 0L, 10L))
    }
    assert(rf.extract(ReferenceRegion("alphabet", 2L, 10L)) === "CTGCGACC")
  }

  sparkTest("construct a set of reference promoters") {
    val g = Gene("theGene",
                 Seq(),
                 true,
                 Iterable(Transcript("transcript1",
                                     Seq(),
                                     "theGene",
                                     true,
                                     Iterable(Exon("exon1",
                                                   "transcript1",
                                                   true,
                                                   ReferenceRegion("alphabet",
                                                                   20L,
                                                                   25L))),
                                     Iterable(),
                                     Iterable())))

    val f = sc.parallelize(Seq(
      "tf1 alphabet 2 4 +",
      "tf2 alphabet 12 14 +",
      "tf3 alphabet 21 23 +"))
    
    val rpArray = ReferencePromoter(sc.parallelize(Seq(g)),
                                    f,
                                    rf,
                                    20,
                                    5,
                                    makeAnnotator(sc)).collect()

    assert(rpArray.length === 1)
    val rp = rpArray.head
    assert(rp.geneId === "theGene")
    assert(rp.pos === ReferenceRegion("alphabet", 0L, 15L))
    assert(rp.tss === ReferencePosition("alphabet", 20L))
    assert(rp.sequence === "AACTGCGACCTATCT")
    assert(rp.tfbs.size === 2)
    assert(rp.tfbs.filter(_.getTf === "tf1").size === 1)
    val tf1 = rp.tfbs.filter(_.getTf === "tf1").head
    assert(tf1.getContig.getContigName === "alphabet")
    assert(tf1.getStart === 2L)
    assert(tf1.getEnd === 5L)
    assert(tf1.getOrientation === Strand.Forward)
    assert(tf1.getSequence === "CTG")
    assert(MathUtils.fpEquals(tf1.getPredictedAffinity, 1.0))
    assert(rp.tfbs.filter(_.getTf === "tf2").size === 1)
    val tf2 = rp.tfbs.filter(_.getTf === "tf2").head
    assert(tf2.getContig.getContigName === "alphabet")
    assert(tf2.getStart === 12L)
    assert(tf2.getEnd === 15L)
    assert(tf2.getOrientation === Strand.Forward)
    assert(tf2.getSequence === "TCT")
    assert(MathUtils.fpEquals(tf2.getPredictedAffinity, 0.5))
  }

  sparkTest("test tfbs load from file") {
    val matchesPath = ClassLoader.getSystemClassLoader
      .getResource("test_matches.txt")
      .getFile
    
    // load text file
    val rdd = sc.textFile(matchesPath)

    val tfbs = ReferencePromoter.extractBindingSites(rdd)
      .map(_._2)
      .collect()

    assert(tfbs.size === 3)
    assert(tfbs.forall(_.getContig.getContigName === "chr1"))
    assert(tfbs.filter(f => f.getStart == 11483L && f.getEnd == 11493L)
      .length === 2)
    assert(tfbs.filter(f => f.getStart == 11482L && f.getEnd == 11494L)
      .length === 1)
    assert(tfbs.map(f => f.getOrientation).distinct.size === 2)
  }
}
