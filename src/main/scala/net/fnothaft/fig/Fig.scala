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
package net.fnothaft.fig

import java.io.File
import net.fnothaft.fig.models._
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.SparkContext._
import org.apache.spark.{ Logging, SparkContext }
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.util.TwoBitFile
import org.bdgenomics.formats.avro._
import org.bdgenomics.utils.cli._
import org.bdgenomics.utils.parquet.io.LocalFileByteAccess
import org.kohsuke.args4j.{ Argument, Option => Args4jOption }
import parquet.avro.AvroReadSupport
import parquet.hadoop.ParquetInputFormat
import parquet.hadoop.util.ContextUtil

object Fig extends BDGCommandCompanion {
  val commandName = "fig"
  val commandDescription = "Find Interesting regulatory Grammar modifying variants"

  def apply(cmdLine: Array[String]) = {
    new Fig(Args4j[FigArgs](cmdLine))
  }
}

class FigArgs extends Args4jBase {
  @Argument(required = true, metaVar = "GENOTYPES", usage = "The genotypes to process.", index = 0)
  var genotypes: String = null

  @Argument(required = true, metaVar = "REFERENCE", usage = "A 2bit file for the reference genome.", index = 1)
  var contigs: String = null

  @Argument(required = true, metaVar = "GENES", usage = "The gene description file.", index = 2)
  var genes: String = null

  @Argument(required = true, metaVar = "MOTIFS", usage = "A file describing TF motifs.", index = 3)
  var motifs: String = null

  @Argument(required = true, metaVar = "SITES", usage = "A feature file describing TFBSs.", index = 4)
  var sites: String = null

  @Argument(required = true, metaVar = "OUTPUT", usage = "Where to write output to.", index = 5)
  var outputPath: String = null

  @Args4jOption(required = false, name = "-before_tss", usage = "The distance to start considering before the TSS. Default is 1500.")
  var startBeforeTss: Int = 1500

  @Args4jOption(required = false, name = "-before_tss", usage = "The distance to stop considering before the TSS. Default is 0.")
  var stopBeforeTss: Int = 0
}

class Fig(protected val args: FigArgs) extends BDGSparkCommand[FigArgs] {
  val companion = Fig

  def run(sc: SparkContext, job: Job) {

    // load two bit file
    val tbf = new TwoBitFile(new LocalFileByteAccess(new File(args.contigs)))

    // load in motifs and create repository
    val motifs = Motif(args.motifs)
    val motifRepository = MotifRepository(sc, motifs)

    // turn these features into structures
    val structures = ReferencePromoter(sc.loadGenes(args.genes),
                                       sc.textFile(args.sites),
                                       tbf,
                                       args.startBeforeTss,
                                       args.stopBeforeTss,
                                       motifRepository)

    // apply variants to promoters
    val variants = VariantPromoter(structures,
                                   sc.loadGenotypes(args.genotypes))

    ???
  }
}
