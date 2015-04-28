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
import org.bdgenomics.adam.models.ReferenceRegion
import org.bdgenomics.formats.avro.Strand

case class LocalMotifRepository(val motifs: Map[String, Seq[Motif]]) extends MotifRepository {
}

// this is a convenience class for testing
class EmptyMotifRepository extends MotifRepository {

  val motifs = Map.empty[String, Seq[Motif]]
  
  override def annotate(tfbs: Iterable[BindingSite],
                        sequence: String,
                        region: ReferenceRegion): Iterable[BindingSite] = {
    Iterable.empty
  }
}
