# [Analyzing Flow Polytopes](#analyzing-flow_polytopes)
Generate and manipulate flow polytopes and their associated graphs. In particular, it generates all possible toric quivers for a given dimension.

## [Getting Started](#getting-started)
This code was written to replicate the results in [\[1\]](#main_paper). Currently, it produces the quivers associated to a given dimension *d*. Further, it can analyze a quiver (directed graph with weights) and give information about the associated polytope. 

The main data structures that are dealt with are Graphs, Quivers, and polytopes at the moment. For graphs and quivers, there are two main representations in python code: 
1. Matrix representation: for a graph <img src="/tex/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92464304999999pt height=22.465723500000017pt/> with edge set <img src="/tex/0e0fff175b21e36dc5c4cae2cb36897c.svg?invert_in_darkmode&sanitize=true" align=middle width=19.477190699999987pt height=22.465723500000017pt/> and vertex set <img src="/tex/b3f35df4c36a139a959bb2514490cd1d.svg?invert_in_darkmode&sanitize=true" align=middle width=19.477190699999987pt height=22.465723500000017pt/>, the matrix associated to <img src="/tex/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92464304999999pt height=22.465723500000017pt/> is the <img src="/tex/62e34e95c7c57ce3b9c407f67922d3e6.svg?invert_in_darkmode&sanitize=true" align=middle width=78.95429354999999pt height=24.65753399999998pt/> matrix <img src="/tex/53d147e7f3fe6e47ee05b88b166bd3f6.svg?invert_in_darkmode&sanitize=true" align=middle width=12.32879834999999pt height=22.465723500000017pt/> with <p align="center"><img src="/tex/af61bd1b82d4f93087079adbbc7d2672.svg?invert_in_darkmode&sanitize=true" align=middle width=273.28447079999995pt height=49.315569599999996pt/></p> Note: for oriented graphs/quivers, we use 1 for the <img src="/tex/3e384b223dce750e6c98aa501355f00b.svg?invert_in_darkmode&sanitize=true" align=middle width=20.679527549999985pt height=21.68300969999999pt/> entry corresponding to the head of an edge and -1 for the tail. 
2. Edge representation: For <img src="/tex/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92464304999999pt height=22.465723500000017pt/> as above, define each edge as a pair <img src="/tex/b31116e7405741aa8eaa5adaa6b31484.svg?invert_in_darkmode&sanitize=true" align=middle width=50.77636904999999pt height=24.65753399999998pt/> of endpoints (ordered by tail/head if it is an oriented graph), and thus we can represent <img src="/tex/5201385589993766eea584cd3aa6fa13.svg?invert_in_darkmode&sanitize=true" align=middle width=12.92464304999999pt height=22.465723500000017pt/> as a list of all such pairs. 

There are 2 major files at the moment:
* **graph\_cal.py** contains routines for generating all of the directed graphs that satisfy conditions to generate a flow polytope in dimension <img src="/tex/2103f85b8b1477f430fc407cad462224.svg?invert_in_darkmode&sanitize=true" align=middle width=8.55596444999999pt height=22.831056599999986pt/>. 
The main routines in this file are: 
    * `Step1(d)`: 
creates all connected undirected graphs up to isomorphism for a given <img src="/tex/2103f85b8b1477f430fc407cad462224.svg?invert_in_darkmode&sanitize=true" align=middle width=8.55596444999999pt height=22.831056599999986pt/> satisfying <img src="/tex/b50f7d37801bbbfa4f855ee8ad3bc566.svg?invert_in_darkmode&sanitize=true" align=middle width=96.43475654999999pt height=24.65753399999998pt/> and <img src="/tex/5bfc9b4f236354b527d35f72b3bedb28.svg?invert_in_darkmode&sanitize=true" align=middle width=137.73828914999999pt height=24.65753399999998pt/> such that the valence of any vertex <img src="/tex/edc38764d2b249d54559b8fd1f3bced4.svg?invert_in_darkmode&sanitize=true" align=middle width=48.12617039999999pt height=22.465723500000017pt/> is <img src="/tex/7e70ada5b76b03379a57d8aca40f296b.svg?invert_in_darkmode&sanitize=true" align=middle width=25.570741349999988pt height=21.18721440000001pt/>. 
    * `Step2(M)`: 
removes the loops from an undirected, connected graph given in matrix form as <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/>.
    * `Step3(M, edges)`: 
takes a graph in matrix form <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> and a subset of the edges in <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> to split by adding a vertex.
    * `Step4(M)`: 
produces all possible oriented quivers <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> for a graph <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/>, such that valence 2 vertices in <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> are sinks and there are no oriented cycles. 
    * `Step5(Ms)`: 
Takes a list <img src="/tex/b1a9dbac8eac63ce661a58d7dfca05f4.svg?invert_in_darkmode&sanitize=true" align=middle width=25.44519449999999pt height=22.465723500000017pt/> of quivers (in matrix form) and returns a list of the ones that are unique(up to directed graph isomorphism)
    * `flow_polytope(Q)`: 
computes the vertices for the convex hull defining the polytope associated to the dual of the quiver <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/>. This algorithm is taken from the procedure outlined in [\[2\]](#neighborly_polytopes) section 3. 
* **quiver\_cal.py** contains routines for studying at the quivers generated in **graph\_cal.py**.  It contains the following routines:
    * `subquivers(M)`: Given a quiver, it generates all its subquivers with the same set of vertices. 
    * `subsets_closed(M)`: 
Give you all the possible subquivers of M such that it is closed under arrows. 
    * `theta(M)`: returns the weights of the vertices.
    * `is_stable(M, subM)`: returns true with a subquiver is stable. 
Note that the input <img src="/tex/fb97d38bcc19230b0acd442e17db879c.svg?invert_in_darkmode&sanitize=true" align=middle width=17.73973739999999pt height=22.465723500000017pt/> is a weighted matrix, so the weights in the arrows


## [Generate Quivers](#generate-quivers)
To create the list of all quivers satisfying requirements for flow polytopes in [\[1\]](#main_paper), can run: 

> `python create_all_possible_graphs.py`

and enter the value for `d` when prompted. 

This generates a set of files in the folder `outputs/d=*/` corresponding to the 4 steps above, together with the list of all quivers, and 2 plots that outline the growth pattern of quivers from initial starting graphs. 

To load these outputs into a list of matrices(matrix representations of all of the different quivers), use the function 
`read_step_file(filename)` from the file `create_all_possible_graphs.py`

## [Flow Polytopes](#flow-polytopes)
At the moment, the code generates the vertices of the dual of the polytope associated to the quiver <img src="/tex/1afcdb0f704394b16fe85fb40c45ca7a.svg?invert_in_darkmode&sanitize=true" align=middle width=12.99542474999999pt height=22.465723500000017pt/> input. Need to change that.

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
