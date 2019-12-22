# Functions in `gc.m2`
* `asList (x)` : turns x into a list 

* `sumList {Axis: None, Row, Col} (x)` : sums the entries in list x 

* `combinations {Replacement, MinSum, MaxSum, Order} (list, num)` : takes combinations of list of length num. 

* `sortedIndices (x)`: returns the order of indices for values in in `sort(x)`

* `isPermutation (x, y)`: answers question: is x a permutation of the rows/cols of y? 

* `unorietnedUniqueUpToPermutation (x)`: returns the list of entries of x that are unique up to permutation, where x is a list of unoriented graphs. 

* `allPossibleBaseGraphsForPair (x)`: the admissable base graphs G with G#0 = x#0 and G#1 = x#1

* `undirectedGraphs (d)`: creates all of the (possibly non-connected) undirected graphs in dimension d. 

* `graphEdges {Oriented: true/false, RavelLoops: true/false} (g)`: returns the edges of graph g as a list of lists. 

* `replaceInList (i, v, l)`: replaces the entry in list at index i with value v

* `graphFromEdges {Oriented: true/false} (E)`: takes list of edges and returns the matrix rep of the graph g

* `isGraphConnected (g)`: answers question: is graph g connected? 

* `edgesOutOfPoint {Oriented} (p, E)`: returns the edges in E that are connected to point p

* `findCycleDFS (startV, visited, E)`: recursive function that, when wrapped, performs DFS search to find cycle in a directed graph given by edge list E

* `existsOrientedCycle (G)`: uses `findCycleDFS` to perform DFS search and answers the question: is there an oriented cycle in the graph G? (G is a directed graph in matrix form)

* `isPathBetween {Oriented, SavePath, EdgesAdded} (p, q, E)`: checks if there exists a path between p and q using edges from E, and optionally returns the path. 

* `isEdgeInCycle (i, E)`: answers question: is there a cycle in E containing edge at index i? 

* `splitLoops (m)`: splits the loops in graph g with matrix rep m. 

* `splitEdges (m, E)`: splits the edges in graph m that are at the indices in list E 

* `Step1 (n)`: returns admissable(connected, unique up to undirected graph isomorphism, and with all edges in a cycle) graphs for quiver generation in dimension n. graphs CAN contain loops. 

* `Step2 (m)`: returns the list of all possible graphs generated from a given output of step1 by splitting subsets of the edges and removing all loops. 

* `Step3 (m)`: add the set of all allowable orientations on an unoriented graph m (note that verticesi of valence 2 must be sinks)

* `Step4 (l)`: returns the values of l (where l is a list of directed graphs) that do not contain any oriented cycles. 

* `Step5 (l)`: returns the values of l (where l is a list of directed graphs) that are unique up to isomorphism. 

* `subquivers (Q)`: returns the subquivers of a quiver Q (given in weighted matrix form)

* `subsetsClosedUnderArrows (Q)`: returns the subquivers of Q that are closed under arrows.

* `theta (Q)`: returns the vertex weights for quiver Q, where the canonical weights on arrows are given.

* `isStable(Q, subQ)`: returns boolean of statement: the arrows corresponding to the list subQ comprise a stable subquiver with respect to the inherited weights from Q

* `unstableSubquivers(Q)`: returns all of the subquivers of Q that are unstable. 

* `isProperSubset(Q1, Q2)`: returns boolean of statement: set(Q1) is a proper subset of set(Q2)

* `maximalUnstableSubquivers(Q)`: returns all of the maximal unstable subquivers of Q

* `isTight(Q)`: returns boolean of statement: Q is a tight quiver
