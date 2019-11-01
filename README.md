# [Analyzing Flow Polytopes](#analyzing-flow_polytopes)
Generate and manipulate flow polytopes and their associated graphs. In particular, it generates all possible toric quivers for a given dimension. Although, currently is overcounting them.

## [Getting Started](#getting-started)
This code was written to replicate the results in [1](#main_paper). Currently, it produces the quivers associated to a given dimension *d*, although not uniquely up to graph isomorphism yet. Further, it can analyze a quiver (directed graph with weights) and give information about the associated polytope. 
There are 2 major files at the moment:
* **graph\_cal.py** contains routines for generating all of the directed graphs that satisfy conditions to generate a flow polytope in dimension *d*
* **quiver\_cal.py** contains routines for studying at the quivers generated in **graph\_cal.py**.  It contains the following routines:

Inpute is a weighted matrix, so the weights in the arrows
subquivers(M): Given a quiver, it generates all its subquivers with the same set of vertices. 
subsets_closed(M): Give you all the possible subquivers of M such that it is closed under arrows. 
theta(M): returns the weights of the vertices.
is_stable(M, subM): returns true with a subquiver is stable. 


## [Generate Quivers](#generate-quivers)

## [Flow Polytopes](#flow-polytopes)

## [References](#references)
<a id='main_paper'>\[1\]
Klaus Altmann, Benjamin Nill, Sabine Schwentner, and Izolda Wiercinska, *Flow polytopes and the graph of reflexive polytopes*, Discrete Mathematics. 309.16(2009), pp 4992-4999. 
[sciencedirect.com/sceince/article/pii/S0012365X09001162](http://www.sciencedirect.com/science/article/pii/S0012365X09001162)</a>
