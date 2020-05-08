# Functions exported from `ToricQuiver.m2'
* `toricQuiver' : construct a Toric Quiver for a range of inputs (connectivity matrix, list of edges)

* `sampleQuiver (d)' : returns a random quiver in dimension d

* `toricQuivers (d)' : returns a list of all d-dimensinoal quivers

* `isTight (Q)': answers question: is the quiver Q tight? (Q given in matrix form)

* `subquivers (Q)' : returns a list of the subquivers of Q, given as sets of vertices 

* `subsetsClosedUnderArrows (Q)': returns the subsets of subquivers(Q) satisfying the closed under arrows condition

* `isStable (Q, sQ)' : answers the question: subQ is a stable subquiver of Q

* `isMaximal (Q, listOfQs)': answers question: Q is maximal with respect to the list of quivers listOfQs

* `isAcyclic (Q)': answers question: Is Q an acyclic quiver?

* `maximalUnstableSubquivers (Q)' : returns the maximal Unstable subqivers of Q

* `theta (Q)' : returns the weights on the vertices for a weighted quiver Q

* `neighborliness (Q)' : returns the neighborliness of quiver Q

* `walls (Q)' : returns the walls for the weight chamber system associated to a weighted quiver Q

* `wallType (Q, Qplus)' : returns the type of wall(wall index set and type (tplus, tminus)) for the wall specified by the partition of Q0 into Qplus and Q0 \ Qplus

* `flowPolytope(Q)' : **under construction** returns the dual polytope associated to Q

* `mergeOnVertex(Q1, v1, Q2, V2)' : **under construction** 

* `mergeOnArrow(Q1, a1, Q2, a2)' : **under construction** 
