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
import net.fnothaft.fig.avro.BindingSite
import org.bdgenomics.adam.models.{ ReferencePosition, ReferenceRegion }
import org.bdgenomics.formats.avro._
import org.bdgenomics.utils.misc.MathUtils
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
                                      Seq(),
                                      new EmptyMotifRepository)
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
                                       refPromoter,
                                       new EmptyMotifRepository)

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
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                          ReferenceRegion("myChrom", 10L, 20L),
                                                          Iterable(Variant.newBuilder()
                                                            .setContig(Contig.newBuilder()
                                                            .setContigName("myChrom")
                                                            .build())
                                                            .setStart(14L)
                                                            .setEnd(15L)
                                                            .setReferenceAllele("A")
                                                            .setAlternateAllele("C")
                                                            .build()))
    assert(haplotype === "ACACCCACAC")
    assert(shifts.isEmpty)
  }

  test("flatten out a haplotype with a deletion") {
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                               ReferenceRegion("myChrom", 10L, 20L),
                                                               Iterable(Variant.newBuilder()
                                                                 .setContig(Contig.newBuilder()
                                                                 .setContigName("myChrom")
                                                                 .build())
                                                                 .setStart(14L)
                                                                 .setEnd(17L)
                                                                 .setReferenceAllele("ACA")
                                                                 .setAlternateAllele("A")
                                                                 .build()))
    assert(haplotype === "ACACACAC")
    assert(shifts.size === 1)
    assert(shifts.head.location === 17L)
    assert(shifts.head.shift === -2)
  }

  test("flatten out a haplotype with an insertion") {
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                               ReferenceRegion("myChrom", 10L, 20L),
                                                               Iterable(Variant.newBuilder()
                                                                 .setContig(Contig.newBuilder()
                                                                 .setContigName("myChrom")
                                                                 .build())
                                                                 .setStart(14L)
                                                                 .setEnd(15L)
                                                                 .setReferenceAllele("A")
                                                                 .setAlternateAllele("ACA")
                                                                 .build()))
    assert(haplotype === "ACACACACACAC")
    assert(shifts.size === 1)
    assert(shifts.head.location === 15L)
    assert(shifts.head.shift === 2)
  }

  test("flatten out a haplotype with a mnp") {
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                               ReferenceRegion("myChrom", 10L, 20L),
                                                               Iterable(Variant.newBuilder()
                                                                 .setContig(Contig.newBuilder()
                                                                 .setContigName("myChrom")
                                                                 .build())
                                                                 .setStart(14L)
                                                                 .setEnd(17L)
                                                                 .setReferenceAllele("ACA")
                                                                 .setAlternateAllele("TCG")
                                                                 .build()))
    assert(haplotype === "ACACTCGCAC")
    assert(shifts.isEmpty)
  }

  test("flatten out a haplotype with two snps") {
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
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
                                                                 .build()))
    assert(haplotype === "ACACTCGCAC")
    assert(shifts.isEmpty)
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
                                        .build()),
                                      new EmptyMotifRepository)

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
                                       refPromoter,
                                       new EmptyMotifRepository)

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
                                       sc.parallelize(Seq(genotype)),
                                       new EmptyMotifRepository)
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
                                       sc.parallelize(genotypes),
                                       new EmptyMotifRepository)
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

  test("annotating a modified site without any tfbs gives an empty iterator") {
    val sites = VariantPromoter.annotateSites("ACACAC",
                                              Iterable.empty,
                                              Iterable.empty,
                                              ReferenceRegion("chrom1", 0L, 6L),
                                              LocalMotifRepository(Map.empty))
    assert(sites.isEmpty)
  }

  test("annotating a modified site whose tfbs can no longer bind gives an empty iterator") {
    val sites = VariantPromoter.annotateSites("ACACAC",
                                              Iterable.empty,
                                              Iterable(BindingSite.newBuilder()
                                                .setTf("tf1")
                                                .setContig(Contig.newBuilder()
                                                .setContigName("chrom1")
                                                .build())
                                                .setStart(1L)
                                                .setEnd(3L)
                                                .setOrientation(Strand.Forward)
                                                .build()),
                                              ReferenceRegion("chrom1", 0L, 6L),
                                              LocalMotifRepository(Map(("tf1" -> Seq(Motif("tf1",
                                                                                           Array(0.0, 1.0, 0.0, 0.0, // C
                                                                                                 0.0, 1.0, 0.0, 0.0) // C
                                                                                         ))))))
    assert(sites.isEmpty)
  }

  test("annotating a modified site whose tfbs can still bind gives us a valid result") {
    val sites = VariantPromoter.annotateSites("ACACAC",
                                              Iterable.empty,
                                              Iterable(BindingSite.newBuilder()
                                                .setTf("tf1")
                                                .setContig(Contig.newBuilder()
                                                .setContigName("chrom1")
                                                .build())
                                                .setStart(1L)
                                                .setEnd(3L)
                                                .setOrientation(Strand.Forward)
                                                .build()),
                                              ReferenceRegion("chrom1", 0L, 6L),
                                              LocalMotifRepository(Map(("tf1" -> Seq(Motif("tf1",
                                                                                           Array(0.0, 1.0, 0.0, 0.0, // C
                                                                                                 0.5, 0.5, 0.0, 0.0) // M
                                                                                         ))))))
    assert(sites.size === 1)
    val site = sites.head
    assert(MathUtils.fpEquals(site.getPredictedAffinity, 0.5))
    assert(site.getTf === "tf1")
    assert(site.getContig.getContigName === "chrom1")
    assert(site.getStart === 1L)
    assert(site.getEnd === 3L)
    assert(site.getOrientation === Strand.Forward)
    assert(site.getSequence === "CA")
  }

  test("we can annotate a modified site that covers a shift") {
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                               ReferenceRegion("chrom1", 10L, 20L),
                                                               Iterable(Variant.newBuilder()
                                                                 .setContig(Contig.newBuilder()
                                                                 .setContigName("chrom1")
                                                                 .build())
                                                                 .setStart(14L)
                                                                 .setEnd(16L)
                                                                 .setReferenceAllele("AC")
                                                                 .setAlternateAllele("A")
                                                                 .build()))
    val sites = VariantPromoter.annotateSites(haplotype,
                                              shifts,
                                              Iterable(BindingSite.newBuilder()
                                                .setTf("tf1")
                                                .setContig(Contig.newBuilder()
                                                .setContigName("chrom1")
                                                .build())
                                                .setStart(14L)
                                                .setEnd(16L)
                                                .setOrientation(Strand.Forward)
                                                .build()),
                                              ReferenceRegion("chrom1", 10L, 20L),
                                              LocalMotifRepository(Map(("tf1" -> Seq(Motif("tf1",
                                                                                           Array(0.75, 0.25, 0.0, 0.0, // M
                                                                                                 0.75, 0.25, 0.0, 0.0) // M
                                                                                         ))))))
    assert(sites.size === 1)
    val site = sites.head
    assert(MathUtils.fpEquals(site.getPredictedAffinity, 0.5625))
    assert(site.getTf === "tf1")
    assert(site.getContig.getContigName === "chrom1")
    assert(site.getStart === 14L)
    assert(site.getEnd === 16L)
    assert(site.getOrientation === Strand.Forward)
    assert(site.getSequence === "AA")
    assert(site.getShift === 0)
  }

  test("we correctly annotate a site that is shifted but not modified") {
    val (haplotype, shifts) = VariantPromoter.extractHaplotype("ACACACACAC",
                                                               ReferenceRegion("chrom1", 10L, 20L),
                                                               Iterable(Variant.newBuilder()
                                                                 .setContig(Contig.newBuilder()
                                                                 .setContigName("chrom1")
                                                                 .build())
                                                                 .setStart(14L)
                                                                 .setEnd(16L)
                                                                 .setReferenceAllele("AC")
                                                                 .setAlternateAllele("A")
                                                                 .build()))
    val sites = VariantPromoter.annotateSites(haplotype,
                                              shifts,
                                              Iterable(BindingSite.newBuilder()
                                                .setTf("tf1")
                                                .setContig(Contig.newBuilder()
                                                .setContigName("chrom1")
                                                .build())
                                                .setStart(16L)
                                                .setEnd(18L)
                                                .setOrientation(Strand.Forward)
                                                .build()),
                                              ReferenceRegion("chrom1", 10L, 20L),
                                              LocalMotifRepository(Map(("tf1" -> Seq(Motif("tf1",
                                                                                           Array(0.75, 0.25, 0.0, 0.0, // M
                                                                                                 0.75, 0.25, 0.0, 0.0) // M
                                                                                         ))))))
    assert(sites.size === 1)
    val site = sites.head
    assert(MathUtils.fpEquals(site.getPredictedAffinity, 0.1875))
    assert(site.getTf === "tf1")
    assert(site.getContig.getContigName === "chrom1")
    assert(site.getStart === 16L)
    assert(site.getEnd === 18L)
    assert(site.getOrientation === Strand.Forward)
    assert(site.getSequence === "AC")
    assert(site.getShift === -1)
  }
}
