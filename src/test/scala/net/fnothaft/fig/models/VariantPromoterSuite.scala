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
import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import org.bdgenomics.formats.avro._
import scala.collection.JavaConversions._
import scala.collection.mutable.Buffer

class VariantPromoterSuite extends FigFunSuite {
  
  def matchesRef(varPromoter: VariantPromoter,
                 refPromoter: ReferencePromoter,
                 sampleName: String,
                 copy: Int,
                 ploidy: Int) {
    assert(varPromoter.promoter === refPromoter)
    assert(varPromoter.sampleId === sampleName)
    assert(varPromoter.copy === copy)
    assert(varPromoter.samplePloidy === ploidy)
    assert(varPromoter.sequence === refPromoter.sequence)
    assert(varPromoter.variants.isEmpty)
  }

  test("if no variants are associated with a promoter, it should match the reference promoter") {
    val refPromoter = ReferencePromoter("myGene",
                                        ReferenceRegion("myChrom", 10L, 20L),
                                        ReferencePosition("myChrom", 25L),
                                        "ACACACACAC",
                                        Iterable())
    val varPromoter = VariantPromoter(refPromoter,
                                      "mySample",
                                      0,
                                      1,
                                      Seq())
    matchesRef(varPromoter,
               refPromoter,
               "mySample",
               0,
               1)
  }

  test("if all the variants for a sample are hom ref, the sample should output two variants that match the reference promoter") {
    val refPromoter = ReferencePromoter("myGene",
                                        ReferenceRegion("myChrom", 10L, 20L),
                                        ReferencePosition("myChrom", 25L),
                                        "ACACACACAC",
                                        Iterable())
    val genotype = Genotype.newBuilder()
      .setVariant(Variant.newBuilder()
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .setContig(Contig.newBuilder()
      .setContigName("myChrom")
      .build())
      .setStart(15L)
      .setEnd(16L)
      .build())
      .setSampleId("mySample")
      .setAlleles(Buffer(GenotypeAllele.Ref, GenotypeAllele.Ref))
      .setIsPhased(true)
      .build()
    val varPromoters = VariantPromoter("myGene",
                                       "mySample",
                                       Seq(genotype),
                                       refPromoter)

    assert(varPromoters.size === 2)
    (0 to 1).foreach(i => {
      val varPromoter = varPromoters.find(_.copy == i)
      assert(varPromoter.isDefined)
      matchesRef(varPromoter.get,
                 refPromoter,
                 "mySample",
                 i,
                 2)
    })    
  }

  test("flatten out a haplotype with a single snp") {
    val (haplotype, _) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                          ReferenceRegion("myChrom", 10L, 20L),
                                                          Iterable(Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(14L)
                                                            .setEnd(15L)
                                                            .setReferenceAllele("A")
                                                            .setAlternateAllele("C")
                                                            .build()),
                                                          Iterable.empty)
    assert(haplotype === "ACACCCACAC")
  }

  test("flatten out a haplotype with a deletion") {
    val (haplotype, _) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                          ReferenceRegion("myChrom", 10L, 20L),
                                                          Iterable(Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(14L)
                                                            .setEnd(17L)
                                                            .setReferenceAllele("ACA")
                                                            .setAlternateAllele("A")
                                                            .build()),
                                                          Iterable.empty)
    assert(haplotype === "ACACACAC")
  }

  test("flatten out a haplotype with an insertion") {
    val (haplotype, _) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                          ReferenceRegion("myChrom", 10L, 20L),
                                                          Iterable(Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(14L)
                                                            .setEnd(15L)
                                                            .setReferenceAllele("A")
                                                            .setAlternateAllele("ACA")
                                                            .build()),
                                                          Iterable.empty)
    assert(haplotype === "ACACACACACAC")
  }

  test("flatten out a haplotype with a mnp") {
    val (haplotype, _) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                          ReferenceRegion("myChrom", 10L, 20L),
                                                          Iterable(Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(14L)
                                                            .setEnd(17L)
                                                            .setReferenceAllele("ACA")
                                                            .setAlternateAllele("TCG")
                                                            .build()),
                                                          Iterable.empty)
    assert(haplotype === "ACACTCGCAC")
  }

  test("flatten out a haplotype with two snps") {
    val (haplotype, _) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                          ReferenceRegion("myChrom", 10L, 20L),
                                                          Iterable(Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(14L)
                                                            .setEnd(15L)
                                                            .setReferenceAllele("A")
                                                            .setAlternateAllele("T")
                                                            .build(), Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(16L)
                                                            .setEnd(17L)
                                                            .setReferenceAllele("A")
                                                            .setAlternateAllele("G")
                                                            .build()),
                                                          Iterable.empty)
    assert(haplotype === "ACACTCGCAC")
  }

  test("if there are variants in a promoter region, they should modify the promoter") {
    val refPromoter = ReferencePromoter("myGene",
                                        ReferenceRegion("myChrom", 10L, 20L),
                                        ReferencePosition("myChrom", 25L),
                                        "ACACACACAC",
                                        Iterable())
    val varPromoter = VariantPromoter(refPromoter,
                                      "mySample",
                                      0,
                                      1,
                                      Seq(Variant.newBuilder()
                                        .setContig(Contig.newBuilder()
                                        .setContigName("myChrom")
                                        .build())
                                        .setStart(14L)
                                        .setEnd(15L)
                                        .setReferenceAllele("A")
                                        .setAlternateAllele("T")
                                        .build()))

    assert(varPromoter.promoter === refPromoter)
    assert(varPromoter.sampleId === "mySample")
    assert(varPromoter.copy === 0)
    assert(varPromoter.samplePloidy === 1)
    assert(varPromoter.sequence === "ACACTCACAC")
    assert(varPromoter.variants.size === 1)
  }

  test("compute variant promoter with a simple genotype") {
    val refPromoter = ReferencePromoter("myGene",
                                        ReferenceRegion("myChrom", 10L, 20L),
                                        ReferencePosition("myChrom", 25L),
                                        "ACACACACAC",
                                        Iterable())
    val genotype = Genotype.newBuilder()
      .setVariant(Variant.newBuilder()
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .setContig(Contig.newBuilder()
      .setContigName("myChrom")
      .build())
      .setStart(14L)
      .setEnd(15L)
      .build())
      .setSampleId("mySample")
      .setAlleles(Buffer(GenotypeAllele.Ref, GenotypeAllele.Alt))
      .setIsPhased(true)
      .build()
    val varPromoters = VariantPromoter("myGene",
                                       "mySample",
                                       Seq(genotype),
                                       refPromoter)

    assert(varPromoters.size === 2)
    val refVarPromoter = varPromoters.find(_.copy == 0)
    assert(refVarPromoter.isDefined)
    matchesRef(refVarPromoter.get,
               refPromoter,
               "mySample",
               0,
               2)

    val altVarPromoter = varPromoters.find(_.copy == 1)
    assert(altVarPromoter.isDefined)
    val varPromoter = altVarPromoter.get
    assert(varPromoter.promoter === refPromoter)
    assert(varPromoter.sampleId === "mySample")
    assert(varPromoter.copy === 1)
    assert(varPromoter.samplePloidy === 2)
    assert(varPromoter.sequence === "ACACCCACAC")
    assert(varPromoter.variants.size === 1)
  }

  sparkTest("compute variant promoter with a simple genotype using RDD methods") {
    val refPromoter = ReferencePromoter("myGene",
                                        ReferenceRegion("myChrom", 10L, 20L),
                                        ReferencePosition("myChrom", 25L),
                                        "ACACACACAC",
                                        Iterable())
    val genotype = Genotype.newBuilder()
      .setVariant(Variant.newBuilder()
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .setContig(Contig.newBuilder()
      .setContigName("myChrom")
      .build())
      .setStart(14L)
      .setEnd(15L)
      .build())
      .setSampleId("mySample")
      .setAlleles(Buffer(GenotypeAllele.Ref, GenotypeAllele.Alt))
      .setIsPhased(true)
      .build()
    val varPromoters = VariantPromoter(sc.parallelize(Seq(refPromoter)),
                                       sc.parallelize(Seq(genotype)))
                                         .collect()

    assert(varPromoters.size === 2)
    val refVarPromoter = varPromoters.find(_.copy == 0)
    assert(refVarPromoter.isDefined)
    matchesRef(refVarPromoter.get,
               refPromoter,
               "mySample",
               0,
               2)

    val altVarPromoter = varPromoters.find(_.copy == 1)
    assert(altVarPromoter.isDefined)
    val varPromoter = altVarPromoter.get
    assert(varPromoter.promoter === refPromoter)
    assert(varPromoter.sampleId === "mySample")
    assert(varPromoter.copy === 1)
    assert(varPromoter.samplePloidy === 2)
    assert(varPromoter.sequence === "ACACCCACAC")
    assert(varPromoter.variants.size === 1)
  }

  sparkTest("compute variant promoter with a simple genotype using RDD methods for multiple samples") {
    val refPromoter = ReferencePromoter("myGene",
                                        ReferenceRegion("myChrom", 10L, 20L),
                                        ReferencePosition("myChrom", 25L),
                                        "ACACACACAC",
                                        Iterable())
    val genotypes = Seq(Genotype.newBuilder()
      .setVariant(Variant.newBuilder()
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .setContig(Contig.newBuilder()
      .setContigName("myChrom")
      .build())
      .setStart(14L)
      .setEnd(15L)
      .build())
      .setSampleId("mySample0")
      .setAlleles(Buffer(GenotypeAllele.Alt, GenotypeAllele.Alt))
      .setIsPhased(true)
      .build(),
                        Genotype.newBuilder()
      .setVariant(Variant.newBuilder()
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .setContig(Contig.newBuilder()
      .setContigName("myChrom")
      .build())
      .setStart(14L)
      .setEnd(15L)
      .build())
      .setSampleId("mySample1")
      .setAlleles(Buffer(GenotypeAllele.Ref, GenotypeAllele.Alt))
      .setIsPhased(true)
      .build(),
                        Genotype.newBuilder()
      .setVariant(Variant.newBuilder()
      .setReferenceAllele("A")
      .setAlternateAllele("C")
      .setContig(Contig.newBuilder()
      .setContigName("myChrom")
      .build())
      .setStart(14L)
      .setEnd(15L)
      .build())
      .setSampleId("mySample2")
      .setAlleles(Buffer(GenotypeAllele.Ref, GenotypeAllele.Ref))
      .setIsPhased(true)
      .build())
    val varPromoters = VariantPromoter(sc.parallelize(Seq(refPromoter)),
                                       sc.parallelize(genotypes))
                                         .collect()

    def checkRef(sampleId: String,
                 copy: Int) {
      val refVarPromoter = varPromoters.filter(_.sampleId == sampleId)
        .find(_.copy == copy)

      assert(refVarPromoter.isDefined)
      matchesRef(refVarPromoter.get,
                 refPromoter,
                 sampleId,
                 copy,
                 2)
    }
    
    def checkAlt(sampleId: String,
                 copy: Int) {
      val altVarPromoter = varPromoters.filter(_.sampleId == sampleId)
        .find(_.copy == copy)

      assert(altVarPromoter.isDefined)
      val varPromoter = altVarPromoter.get
      assert(varPromoter.promoter === refPromoter)
      assert(varPromoter.sampleId === sampleId)
      assert(varPromoter.copy === copy)
      assert(varPromoter.samplePloidy === 2)
      assert(varPromoter.sequence === "ACACCCACAC")
      assert(varPromoter.variants.size === 1)
    }

    assert(varPromoters.size === 6)
    checkAlt("mySample0", 0)
    checkAlt("mySample0", 1)
    checkRef("mySample1", 0)
    checkAlt("mySample1", 1)
    checkRef("mySample2", 0)
    checkRef("mySample2", 1)
  }
}
