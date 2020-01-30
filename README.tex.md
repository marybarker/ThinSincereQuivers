# [Analyzing Flow Polytopes](#analyzing-flow_polytopes)
Generate and manipulate flow polytopes and their associated graphs. In particular, it generates all possible toric quivers for a given dimension.

## [Getting Started](#getting-started)
This code was written to replicate the results in [\[1\]](#main_paper). Currently, it produces the quivers associated to a given dimension *d*. Further, it can analyze a quiver (directed graph with weights) and give information about the associated polytope. 

The main data structures that are dealt with are Graphs, Quivers, and polytopes at the moment. For graphs and quivers, there are two main representations in python code: 
1. Matrix representation: for a graph $G$ with edge set $G_1$ and vertex set $G_0$, the matrix associated to $G$ is the $|G_0|\times|G_1|$ matrix $A$ with $$A_{i,j}=\begin{cases}1&\text{if edge }j\text{ contains vertex }i\\0&\text{ otherwise}\end{cases}$$ Note: for oriented graphs/quivers, we use 1 for the $i, j$ entry corresponding to the head of an edge and -1 for the tail. 
2. Edge representation: For $G$ as above, define each edge as a pair $(v_0, v_1)$ of endpoints (ordered by tail/head if it is an oriented graph), and thus we can represent $G$ as a list of all such pairs. 

There are 2 major files at the moment:
* **graph\_cal.py** contains routines for generating all of the directed graphs that satisfy conditions to generate a flow polytope in dimension $d$. 
The main routines in this file are: 
    * `Step1(d)`: 
creates all connected undirected graphs up to isomorphism for a given $d$ satisfying $|G_0|\le 2d-1$ and $|G_1|=|G_0|+d-1$ such that the valence of any vertex $v\in G_0$ is $\ge 3$. 
    * `Step2(M)`: 
removes the loops from an undirected, connected graph given in matrix form as $M$.
    * `Step3(M, edges)`: 
takes a graph in matrix form $M$ and a subset of the edges in $M$ to split by adding a vertex.
    * `Step4(M)`: 
produces all possible oriented quivers $Q$ for a graph $M$, such that valence 2 vertices in $Q$ are sinks and there are no oriented cycles. 
    * `Step5(Ms)`: 
Takes a list $Ms$ of quivers (in matrix form) and returns a list of the ones that are unique(up to directed graph isomorphism)
    * `flow_polytope(Q)`: 
computes the vertices for the convex hull defining the polytope associated to the dual of the quiver $Q$. This algorithm is taken from the procedure outlined in [\[2\]](#neighborly_polytopes) section 3. 
* **quiver\_cal.py** contains routines for studying at the quivers generated in **graph\_cal.py**.  It contains the following routines:
    * `subquivers(M)`: Given a quiver, it generates all its subquivers with the same set of vertices. 
    * `subsets_closed(M)`: 
Give you all the possible subquivers of M such that it is closed under arrows. 
    * `theta(M)`: returns the weights of the vertices.
    * `is_stable(M, subM)`: returns true with a subquiver is stable. 
Note that the input $M$ is a weighted matrix, so the weights in the arrows


## [Generate Quivers](#generate-quivers)
To create the list of all quivers satisfying requirements for flow polytopes in [\[1\]](#main_paper), can run: 

> `python create_all_possible_graphs.py`

and enter the value for `d` when prompted. 

This generates a set of files in the folder `outputs/d=*/` corresponding to the 4 steps above, together with the list of all quivers, and 2 plots that outline the growth pattern of quivers from initial starting graphs. 

To load these outputs into a list of matrices(matrix representations of all of the different quivers), use the function 
`read_step_file(filename)` from the file `create_all_possible_graphs.py`

## [Flow Polytopes](#flow-polytopes)
At the moment, the code generates the vertices of the dual of the polytope associated to the quiver $Q$ input. Need to change that.

TODO: 
- [x] add cycle detection to M2 version (for step4)
- [ ] change dual to actual polytope
- [ ] visualize polytopes (simple matplotlib implementation)
- [ ] restructure code for clarity(quiver_cal depends on graph_cal)
- [x] add quiver_cal functionality to M2 version

## [References](#references)
<a id='main_paper'>\[1\]
Klaus Altmann, Benjamin Nill, Sabine Schwentner, and Izolda Wiercinska, *Flow polytopes and the graph of reflexive polytopes*, Discrete Mathematics. 309.16(2009), pp 4992-4999. 
[sciencedirect.com/sceince/article/pii/S0012365X09001162](http://www.sciencedirect.com/science/article/pii/S0012365X09001162)</a>

<a id='neighborly_polytopes'>\[2\]
Patricio Gallardo and Daniel Mckenzie, *On the neighborliness of dual flow polytopes of quivers*, 2018, <a href='http://arxiv.org/abs/1811.01993'>arXiv:1811.01993</a>
</a>
