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

    def slice(self, arrows):
        return ToricQuiver(self.connectivity_matrix[:,arrows])


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


def coneIntersection(a, b):
    return np.matrix((np.array(a)[:, None] == np.array(b)).all(-1).any(1))


def intersectionDim(a, b):
    return np.linalg.matrix_rank(np.matrix((a[:, None] == b).all(-1).any(1)))
    #nr, nc = a.shape
    #dtype={'names':['f{}'.format(i) for i in range(nc)],
    #       'formats':nc * [a.dtype]}
    #c = np.intersect1d(a.view(dtype), b.view(dtype))
    #return c.view(a.dtype).reshape(-1, nc)


def coneSystem(Q):
    spanning_trees = allSpanningTrees(Q)
    qcmt = Q.connectivity_matrix.transpose()

    Q_dim = qcmt.shape[0] - qcm.shape[1] + 1
    cone_dim = qcm.shape[1] - 1

    # create list of cones CT for each tree T
    tree_chambers = unique(apply(lambda st: qcmt[st,:], spanning_trees))

    if len(tree_chambers) > 1 then (
        # find all pairs of cones with full-dimensional interesection
        aij = [j \
            for i, tci in enumerate(tree_chambers) \
            for j in range(i+1, len(tree_chambers)) \
            if intersectionDim(tree_chambers[i], 
                               tree_chambers[j]) >= cone_dim \
        ]

        # now add each cone CT to the list of admissable subcones 
        # before adding every possible intersection to the list as well.
        added_to = set(tree_chambers) # keep track of cones that have already been added
        last_list = [(tc, tci) for tci, tc in enumerate(tree_chambers)]
        subsets = list(tree_chambers)

        # generate all possible intersections of cones
        for ctr in range(len(spanning_trees)):
            # stopping condition: if all intersections tried in the latest loop iteration are up lower-dim
            all_empty = True 

            # create a list of intersecting pairs from intersecting treeChamber elements with lastList elements
            current_list = []
            for list_entry in last_list:
                tree_i, i = list_entry
                for j in aij[i].tolist():
                    tree_j = tree_chambers[j]
                    tree_ij = coneIntersection(tree_i, tree_j)
                    if np.linalg.matrix_rank(tree_ij) >= cone_dim:
                        all_empty = False
                        added_to.add(tree_ij)
                        current_list.append((tree_ij, j))
            if all_empty:
               break
            else:
               last_list = list(currentList)
               subsets = subsets + list(zip(*last_list))[0]

        for si in subsets:
            contains_something = False
            for sj in subsets:
                if si != sj and coneIntersection(si, sj).shape[0] >= sj.shape[0]:
                    contains_something = True
                    break
            if not contains_something:
                yield si


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
        f = incInverse(Q.subquiver(t[0]), weight)
        f = incInverse(Q.subquiver(t[0]), weight, False)
        if all(f >= 0):
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


def incInverse(Q, theta, nonneg=True):
    nonzero_flows = [1 if x != 0 else 0 for x in Q.flow]
    nc = Q.connectivity_matrix.shape[1]
    a = np.zeros((nc,nc))
    np.fill_diagonal(a, nonzero_flows)
    a = Q.connectivity_matrix*a

    return gt.intSolve(a, theta, nonneg=nonneg)



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

    subquiver_matrix = (Q.connectivity_matrix[subquiver_vertices,:])[:,SQ]

    for ss in gt.subsetsClosedUnderArrows(subquiver_matrix):
        if sum(weights[list(ss)]) + minimum_weight <= 0:
            return False
    return True


def isTight(Q, flow=None):
    if flow is None:
        flow=Q.flow
    maximal_unstable_subs = maximalUnstableSubquivers(ToricQuiver(Q.connectivity_matrix, flow), return_singletons=True)

    num_arrows = Q.connectivity_matrix.shape[1]
    if num_arrows > 1:
        return all(list(map(lambda x: len(x) != (num_arrows - 1), maximal_unstable_subs["nonsingletons"])))
    else:
        return len(maxUnstSubs["singletons"]) < 1


def makeTight(Q, th):
    if len(th) == len(Q.Q1): # if passing in a flow instead of weights
        potentialF = th
        th = theta(Q, potentialF)
    else:
        potentialF = list(incInverse(Q, th))


    # this function calls itself recursively until a tight quiver is produced
    if isTight(Q, potentialF):
        return ToricQuiver(Q.connectivity_matrix, potentialF);
    else:
        if (len(list(stableTrees(Q, th))) < 1):
            print("Error: provided weight theta is not in C(Q) and so does not admit a tight toric quiver")
            return

        # find a maximal unstable subquiver, called R, of codimension 1 (this exists if not tight)
        dm = np.zeros((len(potentialF), len(potentialF)))
        np.fill_diagonal(dm, potentialF)
        Qcm = matrixFromEdges(Q.Q1)*dm

        max_unstable_subs = maximalUnstableSubquivers(ToricQuiver(Q.connectivity_matrix, potentialF), return_singletons=True)
        R = max_unstable_subs['nonsingletons'][0]
        Rvertices = [x for y in R for x in Q.Q1[y]]
        S = []

        # generate a connected subset of R0 that is closed under arrows and has weight sum <= 0
        if len(R) < 1:
            # this is for the case of quivers with 1 arrow or similar
            Rvertices = max_unstable_subs['singletons'][0]
            S = Rvertices
        else:
            success = False;
            for i in range(1, len(Rvertices)):
                combs = combinations(Rvertices, len(Rvertices)-i)
                for c in combs:
                    if sum([th[x] for x in c]) <= 0:
                        if isClosedUnderArrows(c, Q.connectivity_matrix[:,R]):
                            success = True
                            S = c
                            break
                if success:
                    break

        # alpha is an arrow in Q1 with head in R and tail in Q\R
        alpha = [x for x in range(len(Q.Q1)) if x not in R]
        alpha = alpha[0]
        (aMinus, aPlus) = sorted(Q.Q1[alpha])

        # contract the arrow alpha to create a new quiver
        r1 = range(aMinus)
        r2 = list(range(aMinus+1,len(Q.Q0)))
        r2.remove(aPlus)
        p1 = Q.connectivity_matrix[r1,:]
        p2 = Q.connectivity_matrix[Q.Q1[alpha],:].sum(axis=0).tolist()
        p3 = Q.connectivity_matrix[r2,:]

        new_matrix = np.concatenate((p1,p2,p3))
        new_flow = [f for i, f in enumerate(potentialF) if i != alpha]
        nonempty_edges = np.where(np.absolute(new_matrix).sum(axis=0) > 0)[1]

        new_weight = np.array(np.matmul(new_matrix, np.matrix(potentialF).transpose()).transpose().astype("int32")).ravel()
        new_q = ToricQuiver(new_matrix[:,nonempty_edges])

        return makeTight(new_q, new_flow)


def maxCodimensionUnstable(Q):
    num_arrows = len(Q.Q1)
    maximal_unstables = maximalUnstableSubquivers(Q, True)

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
        for j, subquiver2 in enumerate(nonstable_nonsingletons):
            if not checked[j] and not checked[i]:
                if isProperSubset(subquiver1, subquiver2):
                    checked[i] = True
                    break
    indices = list(filter(lambda x: not checked[x], range(len(checked))))
    with_arrows = [unstable_nonsingletons[y] for y in indices]
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
    checked = [False for s in unstable_nonsingletons]

    for i, subquiver1 in enumerate(unstable_nonsingletons):
        for j, subquiver2 in enumerate(unstable_nonsingletons):
            if not (checked[j] or checked[i] or i == j):
                if isProperSubset(subquiver1, subquiver2):
                    checked[i] = True
                    break

    indices = list(filter(lambda x: not checked[x], range(len(checked))))
    edges = [unstable_nonsingletons[y] for y in indices]
    to_return = {"nonsingletons": edges}

    if return_singletons:
        if len(edges) < 1:
            edges = [[]]

        contained_singletons = [x for y in set(chain(*edges)) for x in Q.Q1[y]]
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


def potentialWalls(Q, theta):
    nv_set = Q.Q0
    num_to_check = int((len(nv_set)+1)/2)

    # create list of possible Qminus
    qminus = [combinations(nv_set, i) for i in range(1, num_to_check)]

    already_met = set()

    walls = []
    for Qm1 in qminus:
        for Qm in Qm1:
            qms = '_'.join([str(x) for x in Qm])
            if not qms in already_met:
                # restrict to only vertices/edgesin Qm
                qm_edges = np.where(np.absolute(Q.connectivity_matrix[Qm,:]).sum(axis=0) == 2)[1]
                Qp = list(set(nv_set).difference(set(Qm)))
                qps = '_'.join([str(x) for x in Qp])
                already_met = already_met.union({qms,qps})

                if gt.isConnected(Qm, [Q.Q1[x] for x in qm_edges]):
                    qp_edges = np.where(np.absolute(Q.connectivity_matrix[Qp,:]).sum(axis=0) == 2)[1]
                
                    if (len(Qp) < 2) or gt.isConnected(Qp, [Q.Q1[x] for x in qp_edges]):
                        walls.append([Qp, wallType(Qp, Q)])
    return walls


def primitiveArrows(Q):
    edges = []
    for i, e in enumerate(Q.Q1):
        other_edges = [Q.Q1[x] for x in range(len(Q.Q1)) if x != i]
        (isCycle, cycle) = gt.existsPathTo(e[0], e[1], other_edges, True, True)
        if not isCycle or len(set(cycle)) < 2:
            edges.append(i)
    return edges


#def referenceThetas(CQ):


def sameChamber(theta1, theta2, CQ):
    """
    this function checks if the weights theta1 and theta2 
    belong to the same chamber in the wall chamber decomposition for Q
    """
    trees_theta1 = stableTrees(Q, theta1)
    trees_theta2 = stableTrees(Q, theta2)
    if len(trees_theta1) < 1:
        if len(trees_theta2) < 1:
            return "cannot be determined. stableTrees are empty"
        else:
            return False
    elif len(trees_theta2) < 1:
         return False
    elif gt.identicalLists(trees_theta1, trees_theta2):
        return true
    else:
        p1 = flowPolytope(Q, weight=theta1)
        p2 = flowPolytope(Q, weight=theta2)
        if p1 == p2:
            return True
        else:
            return "cannot be determined"


def spanningTree(Q, tree_format="edge"):
    return gt.spanningTree(Q.connectivity_matrix, tree_format=tree_format)


def stableTrees(Q, weight):
    for s in allSpanningTrees(Q, tree_format="vertex"):
        ii = incInverse(Q.slice(s[0]), weight)
        if all(incInverse(Q.slice(s[0]), weight) > 0):
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
    if not isinstance(M, np.matrix):
        M = M.connectivity_matrix
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


def wallType(Qplus, Q):
    tplus = np.where(Q.connectivity_matrix[Qplus,:].sum(axis=0) < 0, 1, 0).sum()
    tminus = np.where(Q.connectivity_matrix[Qplus,:].sum(axis=0) > 0, 1, 0).sum()
    return (tplus, tminus)

