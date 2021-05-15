import graph_tools as gt
import numpy as np
from itertools import combinations, chain

edgesFromMatrix = gt.edgesFromMatrix
matrixFromEdges = gt.matrixFromEdges

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
    return gt.allSpanningTrees(Q.connectivity_matrix, tree_format=tree_format)


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
    if isinstance(Q, np.matrix):
        return gt.isClosedUnderArrows(V, Q)
    else:
        return gt.isClosedUnderArrows(V, Q.connectivity_matrix)


def isProperSubset(A, B):
    return set(A) < set(B)


def isSemistable(SQ, Q):
    # get the vertices in the subquiver
    subquiver_vertices = list(set([x for y in SQ for x in Q.Q1[y]]))
    other_vertices = [x for x in Q.Q0 if x not in subquiver_vertices]

    # inherited weights on the subquiver
    #weights = Q.weight[other_vertices]
    weights = np.array([Q.weight[x] for x in subquiver_vertices])

    # negative weights in Q_0 \ subQ_0
    minimum_weight = sum(map(lambda x: x if x <=0 else 0, [0] + weights))

    subquiver_matrix = (Q[SQ])[subquiver_vertices,:]

    for ss in gt.subsetsClosedUnderArrows(subquiver_matrix):
        if sum(weights[list(ss)]) + minimum_weight < 0:
            return False
    return True


def isStable(SQ, Q):
    # get the vertices in the subquiver
    subquiver_vertices = list(set([x for y in SQ for x in Q.Q1[y]]))
    other_vertices = [x for x in Q.Q0 if x not in subquiver_vertices]

    # inherited weights on the subquiver
    weights = np.array([Q.weight[x] for x in subquiver_vertices])

    # negative weights in Q_0 \ subQ_0
    minimum_weight = sum(map(lambda x: x if x <=0 else 0, [0] + weights))

    subquiver_matrix = (Q[SQ])[subquiver_vertices,:]

    for ss in gt.subsetsClosedUnderArrows(subquiver_matrix):
        if sum(weights[list(ss)]) + minimum_weight <= 0:
            return False
    return True


def isTight(Q):
    maximal_unstable_subs = maximalUnstableSubquivers(Q, return_singletons=True)

    num_arrows = Q.connectivity_matrix.shape[1]
    if num_arrows > 1:
        return all(apply(lambda x: len(x) != (num_arrows - 1), maximal_unstable_subs["nonsingletons"]))
    else:
        return len(maxUnstSubs["singletons"]) < 1


#def makeTight(Q, theta):


def maxCodimensionUnstable(Q):
    num_arrows = len(Q.Q1)
    maximal_unstables = maximalUnstableSubquivers(Q, True)
    print("looking at ", maximal_unstables)

    if len(maximal_unstables["singletons"]) > 0:
        return num_arrows
    else:
        max_unst_lens = min(map(lambda x: len(x), maximal_unstables["nonsingletons"]))
        return num_arrows - max_unst_lens


def maximalNonstableSubquivers(Q):
    nonstable_nonsingletons = list(nonstableSubquivers(Q, "list"))

    with_arrows = []
    checked = [False for s in nonstable_nonsingletons]

    for i, subquiver1 in enumerate(nonstable_nonsingletons):
        is_maximal = True
        for j, subquiver2 in enumerate(nonstable_nonsingletons):
            if not checked[j] and not checked[i]:
                if isProperSubset(subquiver1, subquiver2):
                    is_maximal = False
                    checked[i] = True
                    break
        if is_maximal:
            with_arrows.append(subquiver1)
    to_return = {"nonsingletons": with_arrows}

    if return_singletons:
        if len(with_arrows) < 1:
            with_arrows = [[]]

        contained_singletons = [x for y in set(chain(*with_arrows)) for x in Q.Q1[y]]
        nonstable_singletons = set([i for i, x in enumerate(Q.weight) if x <= 0])

        to_return["singletons"] = list(nonstable_singletons.difference(contained_singletons))
    return to_return


def maximalUnstableSubquivers(Q, return_singletons=False):
    unstable_nonsingletons = list(unstableSubquivers(Q, "list"))

    with_arrows = []
    checked = [False for s in unstable_nonsingletons]

    for i, subquiver1 in enumerate(unstable_nonsingletons):
        is_maximal = True
        for j, subquiver2 in enumerate(unstable_nonsingletons):
            if not checked[j] and not checked[i]:
                if isProperSubset(subquiver1, subquiver2):
                    is_maximal = False
                    checked[i] = True
                    break
        if is_maximal:
            with_arrows.append(subquiver1)
    to_return = {"nonsingletons": with_arrows}

    if return_singletons:
        if len(with_arrows) < 1:
            with_arrows = [[]]

        contained_singletons = [x for y in set(chain(*with_arrows)) for x in Q.Q1[y]]
        unstable_singletons = set([i for i, x in enumerate(Q.weight) if x < 0])

        to_return["singletons"] = list(unstable_singletons.difference(contained_singletons))
    return to_return


#def mergeOnArrow(Q1,a1,Q2,a2):


#def mergeOnVertex(Q1,v1,Q2,v2):


def nonstableSubquivers(Q, output_format="subquiver"):
    if output_format=="subquiver":
        for x in subquivers(Q):
            if not isStable(x, Q):
                yield Q.subquiver(x)
    else:
        for x in subquivers(Q):
            if not isStable(x, Q):
                yield x


#def potentialWalls(Q, theta)


def primitiveArrows(Q):
    edges = []
    for i, e in enumerate(Q.Q1):
        other_edges = [Q.Q1[x] for x in range(len(Q.Q1)) if x != i]
        (isCycle, cycle) = gt.existsPathTo(e[0], e[1], other_edges, True, True)
        if not isCycle or len(set(cycle)) < 2:
            edges.append(i)
    return edges


#def referenceThetas(CQ):


#def sameChamber(theta1, theta2, CQ):


def spanningTree(Q, tree_format="edge"):
    return gt.spanningTree(Q.connectivity_matrix, tree_format=tree_format)


def stableTrees(Q, weight):
    Qalt = ToricQuiver(Q, flow=incInverse(Q, weight))
    for s in allSpanningTrees(Q):
        if isStable(s[0], Qalt):
            yield s


def subsetsClosedUnderArrows(Q):
    if isinstance(Q, np.matrix):
        return gt.subsetsClosedUnderArrows(Q)
    else:
        return gt.subsetsClosedUnderArrows(Q.connectivity_matrix)


def subquivers(Q):
    Q1 = len(Q.Q1)
    for i in range(1, Q1):
        for c in combinations(range(Q1), i):
            yield c


def theta(M, F=None):
    if F is None:
        F = np.ones(M.shape[1])
    return np.array(np.matmul(M, np.matrix(F).transpose()).transpose().astype("int32")).ravel()


def threeVertexQuiver(a,b,c, flow="default"):
    edges = [[0,1] for i in range(a)] \
          + [[1,2] for i in range(b)] \
          + [[0,2] for i in range(c)]

    return ToricQuiver(Es, flow=flow)


def unstableSubquivers(Q, output_format="subquiver"):
    if output_format=="subquiver":
        for x in subquivers(Q):
            if not isSemistable(x, Q):
                yield Q.subquiver(x)
    else:
        for x in subquivers(Q):
            if not isSemistable(x, Q):
                yield x


#def wallType(W):
