import numpy as np
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

def edgesFromMatrix(mat):
    return [[r.index(-1), r.index(1)] for r in mat.transpose().tolist()]


def matrixFromEdges(edges, oriented=True):
    nv = len(set(np.array(edges).flatten()))
    tail = -1 if oriented else 1
    return np.matrix([[tail if x == e[0] else 1 if x == e[1] else 0 for x in range(nv)] for e in edges]).transpose()


def theta(M, F=None):
    if F is None:
        F = np.ones(M.shape[1])
    return np.matmul(M, np.matrix(F).transpose()).transpose().astype("int32").tolist()


class ToricQuiver():
    def __init__(self, edges, flow="default"):
        self.Q2 = edges
        self.connectivity_matrix = matrixFromEdges(edges)
        self.Q1 = range(self.connectivity_matrix.shape[0])
        if flow == "default":
            self.flow = list(np.ones(len(edges), dtype="int32"))
        elif isinstance(flow, list) and len(flow) == len(edges):
            self.flow = flow
        else:
            return "Error in ToricQuiver init: provided flow is not compatible with edges"
        self.weight = theta(self.connectivity_matrix, self.flow)

    def __init__(self, connectivity_matrix, flow="default"):
        self.Q1 = range(connectivity_matrix.shape[0])
        self.connectivity_matrix = connectivity_matrix
        self.Q2 = edgesFromMatrix(connectivity_matrix)

        if flow == "default":
            self.flow = list(np.ones(len(edges), dtype="int32"))
        elif isinstance(flow, list) and len(flow) == len(edges):
            self.flow = flow
        else:
            return "Error in ToricQuiver init: provided flow is not compatible with edges"
        self.weight = theta(self.connectivity_matrix, self.flow)


def allSpanningTrees(Q):
    M = Q.connectivity_matrix
    Q0,Q1 = M.shape

    all_edges = edgesFromMatrix(M)
    all_nodes = range(Q0)

    trees = []
    edge_indices = []
    
    d = Q1 - Q0 + 1
    if d > 0 then (
        # try removing every combination of d edges and see if result is a tree
        d_tuples_to_remove = list(combinations(range(Q1), d))
        edges_kept = []
        edges_removed = []

        for d_tuple in d_tuples_to_remove:

            edges_kept = [e for i, e in enumerate(all_edges) if i not in d_tuple]
            edges_removed = [all_edges[d] for d in d_tuple]

            if isConnected(all_nodes, edges_kept) and isAcyclic(edges_kept):
                trees.append([edges_kept, edges_removed])
        return trees
    return [],[]

'''
def basisForFlowPolytope(Q, spanningTree=None):

def bipartiteQuiver(n, m):

def chainQuiver(l):

def coneSystem(Q):

def flowPolytope(Q, weight=None):

def incInverse(Q, theta):

def isAcyclic(Q):

def isClosedUnderArrows(V, Q):

def isSemistable(SQ, Q):

def isStable(SQ, Q):

def isTight(Q):

def makeTight(Q, theta):

def maxCodimensionUnstable(Q):

def maximalNonstableSubquivers(Q):

def maximalUnstableSubquivers(Q):

def mergeOnArrow(Q1,a1,Q2,a2):

def mergeOnVertex(Q1,v1,Q2,v2):

def potentialWalls(Q, theta)

def primitiveArrows(Q):

def referenceThetas(CQ):

def sameChamber(theta1, theta2, CQ):

def stableTrees(Q, weight):

def subquivers(Q):

def theta(Q):

def threeVertexQuiver(a,b,c):

def wallType(W):

quiverConnectivityMatrix
quiverEdges
quiverFlow
quiverWeights

'''
