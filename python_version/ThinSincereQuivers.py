import graph_tools as gt
import numpy as np
from itertools import product, combinations, permutations, combinations_with_replacement

'''
-- Methods/Functions
    "allSpanningTrees",
    "basisForFlowPolytope",
    "bipartiteQuiver",
    "chainQuiver",
    "coneSystem",
    "flowPolytope",
    "incInverse",
    "isAcyclic",
    "isClosedUnderArrows",
    "isSemistable",
    "isStable",
    "isTight",
    "makeTight",
    "maxCodimensionUnstable",
    "maximalNonstableSubquivers",
    "maximalUnstableSubquivers",
    "mergeOnArrow",
    "mergeOnVertex",
    "potentialWalls",
    "primitiveArrows",
    "quiverConnectivityMatrix",
    "quiverEdges",
    "quiverFlow",
    "quiverWeights",
    "referenceThetas",
    "sameChamber",
    "stableTrees",
    "subquivers",
    "theta",
    "threeVertexQuiver",
    "wallType",

    "toricQuiver"
}
'''

matrixFromEdges = gt.matrixFromEdges
edgesFromMatrix = gt.edgesFromMatrix

def theta(M, F=None):
    if F is None:
        F = np.ones(M.shape[1])
    return np.matmul(M, np.matrix(F).transpose()).transpose().astype("int32").tolist()


class ToricQuiver():

    def __init__(self, edges=None, connectivity_matrix=None, flow="default"):
        if edges is not None:
            self.Q1 = edges
            self.connectivity_matrix = matrixFromEdges(edges)
        elif connectivity_matrix is not None:
            self.connectivity_matrix = connectivity_matrix
            self.Q1 = edgesFromMatrix(connectivity_matrix)
        else:
            return "Error in init: must provide either edges or connectivity matrix"
            
        self.Q0 = range(self.connectivity_matrix.shape[0])
        if flow == "default":
            self.flow = list(np.ones(len(edges), dtype="int32"))
        elif isinstance(flow, list) and len(flow) == len(self.Q1):
            self.flow = flow
        else:
            return "Error in ToricQuiver init: provided flow is not compatible with edges"
        self.weight = theta(self.connectivity_matrix, self.flow)

    def __repr__(self):
        return "toric quiver:\nincidence matrix: " \
               +repr(self.connectivity_matrix)+"\nedges: " \
               +repr(self.Q1)+"\nflow: "+repr(self.flow) \
               +"\nweight: "+repr(self.weight)

    def __getitem__(self, i):
        return self.connectivity_matrix[:,i]


def allSpanningTrees(Q, tree_format="edge"):
    M = Q.connectivity_matrix
    Q0,Q1 = M.shape

    all_edges = edgesFromMatrix(M)
    all_nodes = range(Q0)

    trees = []
    edge_indices = []
    
    d = Q1 - Q0 + 1
    if d > 0:
        # try removing every combination of d edges and see if result is a tree
        d_tuples_to_remove = list(combinations(range(Q1), d))
        edges_kept = []
        edges_removed = []

        for d_tuple in d_tuples_to_remove:

            edges_kept = [e for i, e in enumerate(all_edges) if i not in d_tuple]
            edges_removed = [all_edges[d] for d in d_tuple]

            if gt.isConnected(all_nodes, edges_kept) and gt.isAcyclic(matrixFromEdges(edges_kept)):
                if tree_format != "edge":
                    edges_kept = [i for i in range(Q1) if i not in d_tuple]
                    edges_removed = d_tuple
                trees.append([edges_kept, edges_removed])
        return trees
    else:
        if tree_format == "edge":
            return [Q.Q1, []]
        else:
            return [range(Q0),[]]


def basisForFlowPolytope(Q, spanning_tree=None):
    Q1 = len(Q.Q1)
    if spanning_tree is not None:
         spanning_tree = [spanning_tree, [x for x in range(Q1) if x not in spanning_tree]]
    else:
        spanning_tree = spanningTree(Q, tree_format="vertex")

    removed_edges = spanning_tree[1]

    f = []
    for i in removed_edges:
        edge = Q.Q1[i]
        edge_list = spanning_tree[0] + [i]

        cycle = gt.primalUndirectedCycle([Q.Q1[s] for s in spanning_tree[0]] + [edge])

        fi = [0 for k in range(Q1)]
        for j in cycle:
            if j >= 0:
                fi[edge_list[j]] = 1
            else:
                k = -(1 + j)
                fi[edge_list[k]] = -1
        f.append(fi)

    return np.matrix([[f[i][j] for i in range(len(removed_edges))] for j in range(Q1)])


def bipartiteQuiver(n, m, flow="default"):
    return ToricQuiver([[i, j+n] for j in range(m) for i in range(n)], flow=flow)


def chainQuiver(l, flow="default"):
    return ToricQuiver([[i, i+1] for i, li in enumerate(l) for j in range(li)], flow=flow)


#def coneSystem(Q):
'''
def flowPolytope(Q, weight=None, format="simplified_basis"):
    if len(weight) != len(Q.weight):
        print("error: the provided weight is in incorrect dimension")
        return 

    # vertices of flow polytope correspond to regular flows on spanning trees
    all_trees = allSpanningTrees(Q)
    regular_flows = [for t in all_trees]

:= unique for x in allTrees list(
        if all(incInverse(th, Q^x), y -> y >= 0) then (
            incInverse(th, Q^x)
        ) else (continue;)

    -- simplified basis option reduces the dimension to what is strictly necessary.
    -- Recall we can represent the polytope in a lower dimensional subspace of R^Q1, 
    -- since the polytope has dimension |Q1|-|Q0|+1 <= |Q1|
    if instance(opts.Format, String) and (toString(opts.Format) != "SimplifiedBasis") then (
        regularFlows
    ) else (

        -- first generate a basis (ether user-specified or generated from first spanning tree)
        fpb := {};
        if instance(opts.Format, List) then ( -- if Format is a spanning tree
            fpb = basisForFlowPolytope(opts.Format, toricQuiver(Q, incInverse(th, Q)));
        ) else ( -- if Format is the string "SimplifiedBasis" 
            fpb = basisForFlowPolytope(toricQuiver(Q, incInverse(th, Q)));
        );

        -- translate regularFlows to contain origin by subtracting of first of them from all
        kerF := flatten for f in regularFlows list(
            ff := f - regularFlows#0;
            entries transpose solve(fpb, matrix(for fff in ff list ({round fff})))
        );

        -- translate interior point to origin(if interior lattice point exists)
        ip := interiorLatticePoints convexHull transpose matrix kerF;
        if #ip > 0 then (
            ip = first entries transpose first ip;
            for e in kerF list(e - ip)
        ) else (
            kerF
        )
    )
)
flowPolytope ToricQuiver := opts -> Q -> (
    flowPolytope(Q.weights, Q, Format=>opts.Format)
)
'''
#def incInverse(Q, theta):
#def isClosedUnderArrows(V, Q):
#def isSemistable(SQ, Q):
#def isStable(SQ, Q):
#def isTight(Q):
#def makeTight(Q, theta):
#def maxCodimensionUnstable(Q):
#def maximalNonstableSubquivers(Q):
#def maximalUnstableSubquivers(Q):
#def mergeOnArrow(Q1,a1,Q2,a2):
#def mergeOnVertex(Q1,v1,Q2,v2):
#def potentialWalls(Q, theta)
#def primitiveArrows(Q):
#def referenceThetas(CQ):
#def sameChamber(theta1, theta2, CQ):

def spanningTree(Q, tree_format="edge"):
    """ returns the first spanning tree for the ToricQuiver object Q that is found
    """
    M = Q.connectivity_matrix
    Q0,Q1 = M.shape

    all_edges = edgesFromMatrix(M)
    all_nodes = range(Q0)

    edge_indices = []
    
    d = Q1 - Q0 + 1
    if d > 0:
        # try removing every combination of d edges and see if result is a tree
        d_tuples_to_remove = list(combinations(range(Q1), d))
        edges_kept = []
        edges_removed = []

        for d_tuple in d_tuples_to_remove:

            edges_kept = [e for i, e in enumerate(all_edges) if i not in d_tuple]
            edges_removed = [all_edges[d] for d in d_tuple]


            if gt.isConnected(all_nodes, edges_kept) and gt.isAcyclic(matrixFromEdges(edges_kept)):
                if tree_format != "edge":
                    edges_kept = [i for i in range(Q1) if i not in d_tuple]
                    edges_removed = list(d_tuple)
                return edges_kept, edges_removed
    if tree_format == "edge":
        return Q.Q1, []
    else:
        return range(Q0),[]


'''
def stableTrees(Q, weight):

def subquivers(Q):

def threeVertexQuiver(a,b,c):

def wallType(W):
'''
