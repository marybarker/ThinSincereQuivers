import graph_tools as gt
import numpy as np
from itertools import combinations, chain
from scipy.spatial import ConvexHull, HalfspaceIntersection

edgesFromMatrix = gt.edgesFromMatrix
matrixFromEdges = gt.matrixFromEdges

class ToricQuiver():

    def __init__(self, graph_obj, flow="default"):
        if isinstance(graph_obj, list): # if graph input is a list of edges
            self.Q1 = graph_obj
            self.incidence_matrix = matrixFromEdges(graph_obj)
        elif isinstance(graph_obj, np.matrix): # if graph input is an incidence matrix
            self.incidence_matrix = graph_obj
            self.Q1 = edgesFromMatrix(graph_obj)
        else:#isinstance(graph_obj, ToricQuiver):
            self.incidence_matrix = graph_obj.incidence_matrix
            self.Q1 = graph_obj.Q1
            
        self.Q0 = range(self.incidence_matrix.shape[0])
        if flow == "default":
            self.flow = list(np.ones(len(self.Q1), dtype="int32"))
        elif isinstance(flow, list) and len(flow) == len(self.Q1):
            self.flow = flow
        else:
            return "Error in init: provided flow is not compatible with edges"
        self.weight = theta(self.incidence_matrix, self.flow)
        

    def __repr__(self):
        return "toric quiver:\nincidence matrix: " \
               +repr(self.incidence_matrix)+"\nedges: " \
               +repr(self.Q1)+"\nflow: "+repr(self.flow) \
               +"\nweight: "+repr(self.weight)

    def __getitem__(self, i):
        return self.incidence_matrix[:,i]

    def subquiver(self, arrows):
        f = [f if i in arrows else 0 for i, f in enumerate(self.flow)]
        return ToricQuiver(self.incidence_matrix, flow=f)

    def slice(self, arrows):
        return ToricQuiver(self.incidence_matrix[:,arrows])


class Cone():

    def __init__(self, listOfVertices=[]):
        self.dim = 0
        self.eqs = []
        self.hull = None
        self.interior_point = None
        self.tol=0
        self.vertices = []

        if len(listOfVertices) > 0: # if the set of vertices nonempty

            # if there are enough vertices to form an n-simplex: 
            if len(listOfVertices) > len(listOfVertices[0]): 
                self.dim = np.linalg.matrix_rank(np.matrix(listOfVertices))
                self.hull = ConvexHull(listOfVertices)
                # ignore interior points in the convex hull
                self.vertices = [np.array(self.hull.points[x]) for x in self.hull.vertices]
                self.eqs = [(vec[:self.dim], vec[self.dim:]) for vec in self.hull.equations]
                self.interior_point = ((np.matrix(self.vertices)/len(self.vertices)).sum(axis=0)).tolist()[0]
                self.tol = 0.001 * min([np.dot(x - y, x - y)**0.5 for ix, x in enumerate(self.vertices) for y in self.vertices[ix+1:]])


    def contains_point(self, point, tol=0):
        if self.dim < 1:
            return False
        return all([np.dot(n, np.array(point))+o <= tol for (n,o) in self.eqs])


    def joggle_to_interior(self, point, extra_eqs = [], max_its=10000):
        pt = np.array(point)
        all_neg = False
        ctr = 0
        eqs = self.eqs + extra_eqs

        while not all_neg and ctr < max_its:
            ctr +=1; all_neg = True

            for (n,o) in eqs:
                val = np.dot(n, pt) + o
                if val > -1e-10:
                    pt = pt - (self.tol + val)*n
                    all_neg = False
        if all_neg:
            return pt
        return None


    def intersection(self, otherCone):
        if (self.dim != otherCone.dim) or (self.dim < 1) or (otherCone.dim < 1):
            return Cone()
        self_points_in_other = [x for x in self.vertices if otherCone.contains_point(x)]
        other_points_in_self = [x for x in otherCone.vertices if self.contains_point(x)]

        new_pt = None
        if len(self_points_in_other) > len(other_points_in_self) > 1:
            pt = np.matrix(self_points_in_other).sum(axis=0).tolist()[0]
            new_pt = self.joggle_to_interior(pt, extra_eqs = otherCone.eqs)

        elif len(other_points_in_self) > 1:
            pt = np.matrix(other_points_in_self).sum(axis=0).tolist()[0]
            new_pt = otherCone.joggle_to_interior(pt, extra_eqs = otherCone.eqs)

        if new_pt is not None:
            try:
                points = HalfspaceIntersection( \
                        np.matrix(self.hull.equations.tolist() + otherCone.hull.equations.tolist()), \
                        new_pt).intersections
                return Cone(list(points))
            except:
                pass

        return Cone()


    def __repr__(self):
        return ", ".join([str(x) for x in list(self.vertices)])


    def contains_cone(self, other, tol=0):
        return all([self.contains_point(x, tol=tol) for x in other.vertices])


def allSpanningTrees(Q, tree_format="edge"):
    return gt.allSpanningTrees(Q.incidence_matrix, tree_format=tree_format)


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


def coneSystem(Q):
    spanning_trees = allSpanningTrees(Q, tree_format="vertex")
    qcmt = Q.incidence_matrix.transpose()

    Q_dim = qcmt.shape[0] - qcmt.shape[1] + 1
    cone_dim = qcmt.shape[1] - 1

    # create list of cones CT for each tree T
    tree_chambers = list(map(lambda st: qcmt[st[0],:], spanning_trees))
    primitive_arrows = primitiveArrows(Q)

    if len(tree_chambers) > 1:
        # find Vertices in C(Q) that define the subchambers
        vec_to_conebase = -np.matrix(qcmt[primitive_arrows,:]).sum(axis=0).transpose()
        V = set([tuple(r) for mat in tree_chambers for r in mat.tolist()])
        V = [np.array(v) for v in V]

        # find a basis for the lower-dimensional space that is spanned by the vertices in C(Q)
        B = gt.findLowerDimSpace(V)
        A = np.dot(np.linalg.inv(np.dot(B.transpose(),B)),B.transpose())
        lower_dim_conebase = np.dot(A,vec_to_conebase).transpose().tolist()
        lower_dim_trees = [Cone( \
                              np.dot(A,t.transpose()).transpose().tolist()+lower_dim_conebase) \
                          for t in tree_chambers]

        # find all pairs of cones with full-dimensional interesection
        aij = [[i+j+1 \
            for j, tcj in enumerate(lower_dim_trees[i+1:]) \
            if (tci.intersection(tcj).dim >= cone_dim)] \
            for i, tci in enumerate(lower_dim_trees) \
        ]

        # now add each cone CT to the list of admissable subcones 
        # before adding every possible intersection to the list as well.
        #added_to = set([str(t) for t in lower_dim_trees]) 
        last_list = [[tc, tuple([tci])] for tci, tc in enumerate(lower_dim_trees)]
        subsets = list(lower_dim_trees)
        to_drop = set() # short-term deleting to make sure intersections don't get too huge

        if len(aij) > 1:
            subsets = []
            # generate all possible intersections of cones
            for ctr in range(len(spanning_trees)):
                # stopping condition: if all intersections tried in
                # the latest loop iteration are up lower-dim
                all_empty = True 

                # create a list of intersecting pairs from intersecting
                # treeChamber elements with lastList elements
                current_list = []
                for list_entry in last_list:
                    tree_i, ii = list_entry
                    i = ii[0]
                    for j in aij[i]:
                        if j > ii[-1]:
                            tree_j = lower_dim_trees[j]
                            tree_ij = tree_i.intersection(tree_j)
                            if tree_ij.dim >= cone_dim:
                                all_empty = False
                                to_drop.add(ii)
                                current_list.append([tree_ij, ii + tuple([j])])
                if all_empty:
                    break
                last_list = list(current_list)
                subsets = subsets + [x for x in last_list if x[1] not in to_drop]

        # now take only the subsets that form a partition of C(Q), without any overlap
        unique_subs = []
        all_intersections = [set(x[1]) for x in subsets]
        added = set()
        unique_points = set()
        for x in subsets:
            if (x[1] not in added) and not any([set(x[1]) < y for y in all_intersections]):
                added.add(x[1])
                unique_subs.append(np.array(x[0].vertices))
        subsets = unique_subs

        # transform back to higher dimensional space
        subsets = [np.array([y for y in \
                   np.round(np.dot(B, x.transpose())) \
                   if (y.sum() > -1.0e-10)]) \
                   for x in subsets]
        return subsets


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
    nc = Q.incidence_matrix.shape[1]
    a = np.zeros((nc,nc))
    np.fill_diagonal(a, nonzero_flows)
    a = Q.incidence_matrix*a

    return gt.intSolve(a, theta, nonneg=nonneg)



def isClosedUnderArrows(V, Q):
    if isinstance(Q, np.matrix):
        return gt.isClosedUnderArrows(V, Q)
    else:
        return gt.isClosedUnderArrows(V, Q.incidence_matrix)


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

    subquiver_matrix = (Q.incidence_matrix[subquiver_vertices,:])[:,SQ]

    for ss in gt.subsetsClosedUnderArrows(subquiver_matrix):
        if sum(weights[list(ss)]) + minimum_weight <= 0:
            return False
    return True


def isTight(Q, flow=None):
    if flow is None:
        flow=Q.flow
    maximal_unstable_subs = maximalUnstableSubquivers(ToricQuiver(Q.incidence_matrix, flow), return_singletons=True)

    num_arrows = Q.incidence_matrix.shape[1]
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
        return ToricQuiver(Q.incidence_matrix, potentialF);
    else:
        if (len(list(stableTrees(Q, th))) < 1):
            print("Error: provided weight theta is not in C(Q) and so does not admit a tight toric quiver")
            return

        # find a maximal unstable subquiver, called R, of codimension 1 (this exists if not tight)
        dm = np.zeros((len(potentialF), len(potentialF)))
        np.fill_diagonal(dm, potentialF)
        Qcm = matrixFromEdges(Q.Q1)*dm

        max_unstable_subs = maximalUnstableSubquivers(ToricQuiver(Q.incidence_matrix, potentialF), return_singletons=True)
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
                        if isClosedUnderArrows(c, Q.incidence_matrix[:,R]):
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
        p1 = Q.incidence_matrix[r1,:]
        p2 = Q.incidence_matrix[Q.Q1[alpha],:].sum(axis=0).tolist()
        p3 = Q.incidence_matrix[r2,:]

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


def mergeOnArrow(Q1, a1, Q2, a2):
    A1 = Q1.incidence_matrix
    A2 = Q2.incidence_matrix

    # shape of output matrix
    nrow = A1.shape[0] + A2.shape[0] - 2
    ncol = A1.shape[1] + A2.shape[1] - 1

    # get local indices of vertices associated to the edge
    q1_e = Q1.Q1[a1]
    q2_e = Q2.Q1[a2]

    # permute A1 so that column a1 and rows associated to edge a1 are at bottom
    TQ1 = np.column_stack([A1[:,:a1], A1[:,a1+1:], A1[:,a1]])
    TQ1 = np.row_stack([TQ1[x] for x in range(TQ1.shape[0]) if x not in q1_e] + [TQ1[q1_e[0]], TQ1[q1_e[1]]])

    # permute A2 so that column a2 and rows associated to edge a2 are at top
    TQ2 = np.column_stack([A2[:,:a2], A2[:,a2+1:]])
    TQ2 = np.row_stack([TQ2[q2_e[0]], TQ2[q2_e[1]]] + [TQ2[x] for x in range(TQ2.shape[0]) if x not in q2_e])

    A = np.zeros((nrow, ncol))
    A[:A1.shape[0], :A1.shape[1]] = TQ1
    A[A1.shape[0]-2:,A1.shape[1]:] = TQ2
    return A


def mergeOnVertex(Q1, v1, Q2, v2):
    A1 = Q1.incidence_matrix
    A2 = Q2.incidence_matrix

    # shape of output matrix
    nrow = A1.shape[0] + A2.shape[0] - 1
    ncol = A1.shape[1] + A2.shape[1]

    # permute input matrices so row a1 is at bottom of A1, and 
    # row a2 is at top of A1.easier to merge that way
    TQ1 = np.row_stack([A1[:v1], A1[v1+1:], A1[v1]])
    TQ2 = np.row_stack([A2[v2], A2[:v2], A2[v2+1:]])

    # create dummy matrix to hold output
    A = np.zeros((nrow, ncol))
    A[:A1.shape[0], :A1.shape[1]] = TQ1
    A[A1.shape[0]-1:, A1.shape[1]:] = TQ2
    
    return A


def nonstableSubquivers(Q, output_format="subquiver"):
    if output_format=="subquiver":
        for x in subquivers(Q):
            if not isStable(x, Q):
                yield Q.subquiver(x)
    else:
        for x in subquivers(Q):
            if not isStable(x, Q):
                yield x


def potentialWalls(Q):
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
                qm_edges = np.where(np.absolute(Q.incidence_matrix[Qm,:]).sum(axis=0) == 2)[1]
                Qp = list(set(nv_set).difference(set(Qm)))
                qps = '_'.join([str(x) for x in Qp])
                already_met = already_met.union({qms,qps})

                if gt.isConnected(Qm, [Q.Q1[x] for x in qm_edges]):
                    qp_edges = np.where(np.absolute(Q.incidence_matrix[Qp,:]).sum(axis=0) == 2)[1]
                
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


def referenceThetas(CQ):
    """
    input: CQ is a list of matrices m, where the columns of m correspond to vertices 
            defining the cone corresponding to that sub-chamber
    output: list of interior lattice points, one for each chamber
    
    """
    rts = []
    for rays in CQ:
        avg = np.array(rays.sum(axis=1)).reshape(-1)
        rts.append(avg/max(1, np.gcd.reduce(avg)))
    return rts


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
    return gt.spanningTree(Q.incidence_matrix, tree_format=tree_format)


def stableTrees(Q, weight):
    for s in allSpanningTrees(Q, tree_format="vertex"):
        ii = incInverse(Q.slice(s[0]), weight)
        if all(incInverse(Q.slice(s[0]), weight) > 0):
            yield s


def subsetsClosedUnderArrows(Q):
    if isinstance(Q, np.matrix):
        return gt.subsetsClosedUnderArrows(Q)
    else:
        return gt.subsetsClosedUnderArrows(Q.incidence_matrix)


def subquivers(Q):
    Q1 = len(Q.Q1)
    for i in range(1, Q1):
        for c in combinations(range(Q1), i):
            yield c


def theta(M, F=None):
    if not isinstance(M, np.matrix):
        M = M.incidence_matrix
    if F is None:
        F = np.ones(M.shape[1])
    return np.array(np.matmul(M, np.matrix(F).transpose()).transpose().astype("int32")).ravel()


def threeVertexQuiver(a,b,c, flow="default"):
    edges = [[0,1] for i in range(a)] \
          + [[1,2] for i in range(b)] \
          + [[0,2] for i in range(c)]

    return ToricQuiver(edges, flow=flow)


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
    tplus = np.where(Q.incidence_matrix[Qplus,:].sum(axis=0) < 0, 1, 0).sum()
    tminus = np.where(Q.incidence_matrix[Qplus,:].sum(axis=0) > 0, 1, 0).sum()
    return (tplus, tminus)

