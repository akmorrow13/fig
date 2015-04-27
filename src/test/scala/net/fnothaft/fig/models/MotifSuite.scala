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

import org.bdgenomics.utils.misc.MathUtils
import org.scalatest.FunSuite

class MotifSuite extends FunSuite {

  test("motif PWM length must be a multiple of 4") {
    intercept[IllegalArgumentException] {
      Motif("myMotif", Array(0.1, 0.2, 0.7))
    }
  }

  test("values in PWM must be probabilities") {
    intercept[IllegalArgumentException] {
      Motif("myMotif", Array(0.1, 0.2, 1.1, -0.4))
    }
  }

  test("values in PWM must sum to 1.0 per position") {
    intercept[IllegalArgumentException] {
      Motif("myMotif", Array(0.1, 0.2, 0.6, 0.0))
    }
  }

  test("score a string versus a motif") {
    val m = Motif("myMotif", Array(0.5, 0.25, 0.15, 0.1, 0.1, 0.15, 0.25, 0.5))

    assert(MathUtils.fpEquals(m.sequenceProbability("AT"), 0.25))
    assert(MathUtils.fpEquals(m.sequenceProbability("TA"), 0.01))
  }

  test("load motifs from a file") {
    val ap1Path = ClassLoader.getSystemClassLoader
      .getResource("AP1_motifs.txt")
      .getFile

    val motifs = Motif(ap1Path)

    assert(motifs.length === 2)
    assert(motifs.filter(_.label == "AP1_disc1").length === 1)
    assert(motifs.filter(_.label == "AP1_disc1")
      .head
      .pwm
      .length === 40)
    assert(motifs.filter(_.label == "AP1_disc2").length === 1)
    assert(motifs.filter(_.label == "AP1_disc2")
      .head
      .pwm
      .length === 36)
  }
}
