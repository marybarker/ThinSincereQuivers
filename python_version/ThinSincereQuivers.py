import graph_tools as gt
import numpy as np
from itertools import combinations

matrixFromEdges = gt.matrixFromEdges
edgesFromMatrix = gt.edgesFromMatrix

class ToricQuiver():

    def __init__(self, graph_obj, flow="default"):
        if isinstance(graph_obj, list): # if graph input is a list of edges
            self.Q1 = graph_obj
            self.connectivity_matrix = matrixFromEdges(graph_obj)
        elif isinstance(graph_obj, np.matrix): # if graph input is an incidence matrix
            self.connectivity_matrix = graph_obj
            self.Q1 = edgesFromMatrix(graph_obj)
        else:#isinstance(graph_obj, ToricQuiver):
            self.connectivity_matrix = graph_obj.connectivity_matrix
            self.Q1 = graph_obj.Q1
            
        self.Q0 = range(self.connectivity_matrix.shape[0])
        if flow == "default":
            self.flow = list(np.ones(len(self.Q1), dtype="int32"))
        elif isinstance(flow, list) and len(flow) == len(self.Q1):
            self.flow = flow
        else:
            return "Error in init: provided flow is not compatible with edges"
        self.weight = theta(self.connectivity_matrix, self.flow)
        

    def __repr__(self):
        return "toric quiver:\nincidence matrix: " \
               +repr(self.connectivity_matrix)+"\nedges: " \
               +repr(self.Q1)+"\nflow: "+repr(self.flow) \
               +"\nweight: "+repr(self.weight)

    def __getitem__(self, i):
        return self.connectivity_matrix[:,i]

    def subquiver(self, arrows):
        f = [f if i in arrows else 0 for i, f in enumerate(self.flow)]
        return ToricQuiver(self.connectivity_matrix, flow=f)


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


def flowPolytope(Q, weight=None, polytope_format="simplified_basis"):
    if weight is not None:
        if len(weight) != len(Q.weight):
            print("error: the provided weight is in incorrect dimension")
            return 
    else:
        weight=Q.weight

    # vertices of flow polytope correspond to regular flows on spanning trees
    all_trees = allSpanningTrees(Q, tree_format="vertex")
    regular_flows = []
    for t in all_trees:
        f = np.round(np.array(incInverse(Q.subquiver(t[0]), weight)))
        if all(f >= -0):
            regular_flows.append(f)

    # simplified basis option reduces the dimension to what is strictly necessary.
    # Recall we can represent the polytope in a lower dimensional subspace of R^Q1, 
    # since the polytope has dimension |Q1|-|Q0|+1 <= |Q1|
    if isinstance(polytope_format, str) and (str(polytope_format) != "simplified_basis"):
        return np.array(regular_flows, dtype='int32')
    else:
        # first generate a basis (ether user-specified or generated from first spanning tree)
        fpb = []
        if isinstance(polytope_format, list): # if Format is a spanning tree
            fpb = np.linalg.pinv(basisForFlowPolytope(ToricQuiver(Q, flow=list(incInverse(Q, weight))), polytope_format))
        else: # if Format is the string "SimplifiedBasis" 
            fpb = np.linalg.pinv(basisForFlowPolytope(ToricQuiver(Q, flow=list(incInverse(Q, weight)))))

        # translate regularFlows to contain origin by subtracting of first of them from all
        kerF = np.matrix([x - regular_flows[0] for x in regular_flows]).transpose()
        return np.array(np.round(fpb*kerF), dtype='int32')


def incInverse(Q, theta):
    nonzero_flows = [1 if x != 0 else 0 for x in Q.flow]
    nc = Q.connectivity_matrix.shape[1]
    a = np.zeros((nc,nc))
    np.fill_diagonal(a, nonzero_flows)
    return np.array((np.linalg.pinv(Q.connectivity_matrix*a)*(np.matrix(theta).transpose()))).ravel()


def isClosedUnderArrows(V, Q):
    return np.all(Q.connectivity_matrix[V,:].sum(axis=0) >= 0)


def isSemistable(SQ, Q):
    # get the vertices in the subquiver
    subquiver_vertices = list(set([x for y in SQ for x in Q.Q1[y]]))
    other_vertices = [x for x in Q.Q0 if x not in subquiver_vertices]

    # inherited weights on the subquiver
    weights = Q.weights[other_vertices]

    # negative weights in Q_0 \ subQ_0
    minimum_weight = sum(apply(lambda x: x if x <=0 else 0, [0] + weights))

    subquiver_matrix = (Q[SQ])[subquiver_vertices,:]

    for ss in subsetsClosedUnderArrows(subquiver_matrix):
        if sum(weights[ss]) + minimum_weight < 0:
            return False
    return True


def isStable(SQ, Q):
    # get the vertices in the subquiver
    subquiver_vertices = list(set([x for y in SQ for x in Q.Q1[y]]))
    other_vertices = [x for x in Q.Q0 if x not in subquiver_vertices]

    # inherited weights on the subquiver
    weights = Q.weights[other_vertices]

    # negative weights in Q_0 \ subQ_0
    minimum_weight = sum(apply(lambda x: x if x <=0 else 0, [0] + weights))

    subquiver_matrix = (Q[SQ])[subquiver_vertices,:]

    for ss in subsetsClosedUnderArrows(subquiver_matrix):
        if sum(weights[ss]) + minimum_weight <= 0:
            return False
    return True

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


#def stableTrees(Q, weight):
#def subquivers(Q):

def theta(M, F=None):
    if F is None:
        F = np.ones(M.shape[1])
    return np.matmul(M, np.matrix(F).transpose()).transpose().astype("int32").tolist()


#def threeVertexQuiver(a,b,c):
#def wallType(W):
