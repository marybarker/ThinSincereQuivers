# Functions exported from `ToricQuiver.m2'
* `allSpanningTrees(Q)' : returns all of the spanning trees for a quiver

* `basisForFlowPolytope(Q)' : returns a list of vectors that form a basis for the subspace of $R^{Q_0}$ that contains the flow polytope associated to $Q$.

* `bipartiteQuiver(n, m)' : returns the toric quiver on an underlying bipartite graph with $n$ sources and $m$ sinks.o

* `chainQuiver(L)' : returns the quiver on a graph with the form of a chain of length $\#L$, where links between pairs of vertices given by values in the list $L$. 

* `coneSystem(Q)' : returns the cones that partition the space of weights for the quiver $Q$

* `flowPolytope(Q)' : returns the dual polytope associated to Q if Q is stable. 

* `incInverse(th, Q)' : returns a preimage of the weights th under the inc map associated to Q.

* `isAcyclic(Q)' : answers the question: is Q acyclic? 

* `isClosedUnderArrows (V, Q)': answers question: for V a subset of the vertices of Q, is V closed under arrows? 

* `isSemistable (sQ, Q)' : answers the question: is sQ a semi-stable subquiver of Q? 

* `isStable (sQ, Q)' : answers the question: is sQ a stable subquiver of Q? 

* `isTight (Q)' : answers question: is the quiver Q tight? (Q given in matrix form)

* `makeTight(W, Q) : make the quiver Q tight with respect to the weight W.

* `maximalUnstableSubquivers (Q)' : returns the maximal Unstable subqivers of Q

* `mergeOnArrow(Q1, a1, Q2, a2)' : Creates a new quiver by identifying arrow a1 from quiver Q1 with arrow a2 of quiver Q2. 

* `mergeOnVertex(Q1, v1, Q2, v2)' : Creates a new quiver by identifying vertex v1 from quiver Q1 with vertex v2 of quiver Q2. 

* `neighborliness (Q)' : returns the neighborliness of quiver Q

* `potentialWalls (Q)' : returns the walls for the weight chamber system associated to a weighted quiver Q

* `primitiveArrows (Q)' : return a list of the indices of edges in Q that are primitive

* `referenceThetas (CS)' : returns a candidate weight in each cone of the list of cones CS

* `sameChamber (TH1, TH2, Q)': answers question: is TH1 in the same chamber as TH2 for the chamber decomposition of weights associated to Q

* `stableTrees (TH, Q)' : returns a list of the spanning trees of Q that are TH-stable

* `subquivers (Q)' : returns a list of the subquivers of Q, given as sets of vertices 

* `theta (Q)' : returns the weights on the vertices for a weighted quiver Q

* `threeVertexQuiver(L)' : return the triangular quiver with three vertices and the number of edges between each pair of vertices specified by values of the list $L$.

* `toricQuiver' : construct a Toric Quiver for a range of inputs (connectivity matrix, list of edges, etc.)

* `wallType (Qplus, Q)' : returns the type of wall(wall index set and type (tplus, tminus)) for the wall specified by the partition of Q0 into Qplus and Q0 \ Qplus

