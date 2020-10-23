# Functions exported from `ToricQuiver.m2'
* `bipartiteQuiver(n, m)' : returns the toric quiver on an underlying bipartite graph with $n$ sources and $m$ sinks.o

* `chainQuiver(L)' : returns the quiver on a graph with the form of a chain of length $\#L$, where links between pairs of vertices given by values in the list $L$. 

* `flowPolytope(Q)' : returns the dual polytope associated to Q if Q is stable. 

* `incInverse(Q, th)' : returns a preimage of the weights th under the inc map associated to Q.

* `isAcyclic(Q)' : answers the question: is Q acyclic? 

* `isClosedUnderArrows (Q, V)': answers question: for V a subset of the vertices of Q, is V closed under arrows? 

* `isSemistable (Q, sQ)' : answers the question: is sQ a semi-stable subquiver of Q? 

* `isStable (Q, sQ)' : answers the question: is sQ a stable subquiver of Q? 

* `isTight (Q)' : answers question: is the quiver Q tight? (Q given in matrix form)

* `makeTight(Q, W) : make the quiver Q tight with respect to the weight W.

* `maximalUnstableSubquivers (Q)' : returns the maximal Unstable subqivers of Q

* `mergeOnArrow(Q1, a1, Q2, a2)' : Creates a new quiver by identifying arrow a1 from quiver Q1 with arrow a2 of quiver Q2. 

* `mergeOnVertex(Q1, v1, Q2, v2)' : Creates a new quiver by identifying vertex v1 from quiver Q1 with vertex v2 of quiver Q2. 

* `neighborliness (Q)' : returns the neighborliness of quiver Q

* `subquivers (Q)' : returns a list of the subquivers of Q, given as sets of vertices 

* `theta (Q)' : returns the weights on the vertices for a weighted quiver Q

* `threeVertexQuiver(L)' : return the triangular quiver with three vertices and the number of edges between each pair of vertices specified by values of the list $L$.

* `toricQuiver' : construct a Toric Quiver for a range of inputs (connectivity matrix, list of edges, etc.)

* `toricQuivers (d)' : returns a list of all d-dimensinoal quivers

* `walls (Q)' : returns the walls for the weight chamber system associated to a weighted quiver Q

* `wallType (Q, Qplus)' : returns the type of wall(wall index set and type (tplus, tminus)) for the wall specified by the partition of Q0 into Qplus and Q0 \ Qplus


