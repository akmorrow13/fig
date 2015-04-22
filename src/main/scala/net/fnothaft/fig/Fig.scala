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
import org.apache.hadoop.mapreduce.Job
import org.apache.spark.SparkContext._
import org.apache.spark.{ Logging, SparkContext }
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.adam.rdd.ADAMContext._
import org.bdgenomics.adam.rdd.BroadcastRegionJoin
import org.bdgenomics.formats.avro._
import org.bdgenomics.utils.cli._
import org.bdgenomics.utils.parquet.rdd.BDGParquetContext._
import org.kitesdk.data.{ DatasetDescriptor, Datasets, Formats, View }
import org.kitesdk.data.mapreduce.DatasetKeyOutputFormat
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

  @Argument(required = true, metaVar = "CONTIGS", usage = "The reference genome.", index = 1)
  var contigs: String = null

  @Argument(required = true, metaVar = "GENES", usage = "The gene description file.", index = 2)
  var genes: String = null

  @Args4jOption(required = false, name = "-before_tss", usage = "The distance to start considering before the TSS. Default is 1500.")
  var startBeforeTss: Int = 1500

  @Args4jOption(required = false, name = "-before_tss", usage = "The distance to stop considering before the TSS. Default is 0.")
  var stopBeforeTss: Int = 0
}

class Fig(protected val args: FigArgs) extends BDGSparkCommand[FigArgs] {
  val companion = Fig

  def run(sc: SparkContext, job: Job) {
    ???
  }
}
