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
import scala.annotation.tailrec
import scala.io.Source

object Motif {

  def apply(filepath: String): Seq[Motif] = {
    // load in file
    val motifSource = scala.io.Source.fromFile(filepath)
    val motifs = motifSource.mkString
    motifSource.close()

    // motif descriptions start with ">"
    motifs.split('>')
      .map(motifBlock => {
        // split motif block at each newline
        val motifLines = motifBlock.split('\n')

        // motif name is in first line
        val motifName = motifLines.head
          .takeWhile(_ != ' ')

        // probability lines are formatted "%c %1.6f %1.6f %1.6f %1.6f"
        val probabilities = motifLines.drop(1).flatMap(l => {
          val a = l.drop(2).take(8).toDouble
          val c = l.drop(11).take(8).toDouble
          val g = l.drop(20).take(8).toDouble
          val t = l.drop(29).take(8).toDouble
          Array(a, c, g, t)
        }).toArray

        Motif(motifName, probabilities)
      })
    .filter(m => m.label.length > 0 && m.pwm.length > 0)
  }
}

case class Motif(label: String,
                 pwm: Array[Double]) {
  require(pwm.length % 4 == 0, "Position weight matrix length must be a multiple of 4.")
  val length = pwm.length / 4

  // each value must be between 0 and 1
  pwm.foreach(p => require(p >= 0.0 && p <= 1.0, "P = %f must be between [0, 1].".format(p)))

  // check that each row sums to 1.0
  (0 until length).foreach(i => {
    val idx = i * 4
    val p = pwm(idx) + pwm(idx + 1) + pwm(idx + 2) + pwm(idx + 3)
    require(MathUtils.fpEquals(p, 1.0),
            "Probability (p = %f) of row %d was not equal to 1.".format(p, i))
  })

  private def baseToIdx(c: Char): Int = c match {
    case 'A' => 0
    case 'C' => 1
    case 'G' => 2
    case 'T' => 3
    case _ => throw new IllegalStateException("Invalid character %s.".format(c))
  }

  def sequenceProbability(sequence: String): Double = {
    @tailrec def score(seq: Iterator[Char],
                       pos: Int,
                       p: Double): Double = {
      // continue until we run out of bases
      // if p = 0, we can terminate early
      if (!seq.hasNext || p == 0.0) {
        p
      } else {
        score(seq, pos + 1, p * pwm(pos * 4 + baseToIdx(seq.next)))
      }
    }

    score(sequence.toIterator, 0, 1.0)
  }
}
