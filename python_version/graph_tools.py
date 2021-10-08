import numpy as np
import cvxpy as cp
from itertools import combinations
from scipy.linalg import qr as QR
#from pyhull.convex_hull import ConvexHull


def allSpanningTrees(M, tree_format="edge"):
    Q0,Q1 = M.shape

    all_edges = edgesFromMatrix(M)
    all_nodes = range(Q0)

    trees = []
    edge_indices = []
    
    d = Q1 - Q0 + 1
    if d > 0:
        # try removing every combination of d edges and see if result is a tree
        d_tuples_to_remove = combinations(range(Q1), d)
        edges_kept = []
        edges_removed = []

        for d_tuple in d_tuples_to_remove:

            edges_kept = [e for i, e in enumerate(all_edges) if i not in d_tuple]
            edges_removed = [all_edges[d] for d in d_tuple]

            if isConnected(all_nodes, edges_kept) and isAcyclic(matrixFromEdges(edges_kept)):
                if tree_format != "edge":
                    edges_kept = [i for i in range(Q1) if i not in d_tuple]
                    edges_removed = d_tuple
                trees.append([edges_kept, edges_removed])
        return trees
    else:
        if tree_format == "edge":
            return [all_edges, []]
        else:
            return [range(Q1),[]]


def edgesFromMatrix(mat):
    return [[r.index(-1), r.index(1)] for r in np.matrix(mat).transpose().tolist()]


def edgesOutOf(p, edge_list, oriented=False):
    """ get all edges that are connected to the vertex p
        inputs:
            - p(int): vertex to get the edges emanating from
            - edge_list(list of tuples): all edges in the graph
            - oriented(bool): whether or not the ordering (a, b) of tuples in edge_list matters. 
              * if True: then returns a list of all tuples of the form (p, v) (where v is any vertex) that occurs in edge_list
              * if False: then returns a list of all tuples of the form (p, v) or (v, p) where... 
        outputs:
           sublist of edge_list of the edges that are either connected to or else coming out of p
    """
    if oriented:
        return [(i, e) for i, e in enumerate(edge_list) if p == e[0]]
    return [(i, e) for i, e in enumerate(edge_list) if p in e]


# check if there exists a cycle in a (possibly unconnected)
# oriented graph, passed in matrix form. 
def existsOrientedCycle(M):
    E = edgesFromMatrix(M)
    v = M.shape[0]

    for first_v in range(v):
        visited = [0 if x != first_v else 1 for x in range(v)]
        result = findCycleDFS(first_v, visited, E)
        if result:
            return True
    return False


def existsPathTo(p, q, edge_list, oriented=False, returnPath=False, savedPath=[]):
    """ checks if there exists a path from vertex p to vertex q that 
        can be obtained by stringing together edges from edge_list. 
        optionally saves and returns the path of edges as well.
        inputs:
            - p(int): starting vertex to find path from
            - q(int): ending vertex to find path to
            - edge_list(list): list of tuples (v1, v2) corresponding to edges containing vertices
            - oriented(bool): T/F of statement "ordering of vertices as edge endpoints matters"
            - returnPath(bool): whether or not to return the path between p and q
            - savedPath(list): internal storage for constructing the path for returnPath=True
        outputs:
            - is_path(bool): True/False of statement "there exists a path from p to q
            - currentPath(optional, list): edges of path between p and q
    """
    is_path = False
    currentPath = []

    for edge in edgesOutOf(p, edge_list, oriented=oriented):
        i, e = edge
        # extract endpoints of edge
        v = e[0]
        if p == e[0]:
            v = e[1]

        if returnPath:
            currentPath = savedPath + [(p, v)]

        # check if there's an edge between p & q
        if q == v:
             is_path = True
             break
        else:
            path = []
            remaining_edges = edge_list[:i] + edge_list[i + 1:]
            if returnPath:
                success, path = existsPathTo(v, q, remaining_edges, oriented, returnPath, savedPath)
            else:
                success = existsPathTo(v, q, remaining_edges, oriented, False)
            if success:
                is_path = True
                currentPath += path
                break
    if returnPath:
        return is_path, currentPath
    return is_path



#DFS search to find cycle in directed graph:
def findCycleDFS(start_v, visited, E):
    ret_val = False
    edges_out = edgesOutOf(start_v, E, oriented=True)

    for edge in edges_out:
        current_visited = list(visited)
        edge_verts = edge[1]
        end_v = edge_verts[1]
        if visited[end_v] == 1:
            return True

        current_visited = visited[:end_v] + [1]
        if end_v < len(visited) - 1:
            current_visited = current_visited + visited[end_v+1:]
        ret_val = findCycleDFS(end_v, current_visited, E)
        if ret_val:
            return True
    return ret_val


def findLowerDimBasis(points):
    # find the (hopefully lower dimensional) subspace 
    # that contains the set of points, which should be input
    # in the format of a list of lists

    if len(points) < len(points[0]):
        return np.matrix(points).transpose()

    n = np.linalg.matrix_rank(np.matrix(points))
    q,r = QR(np.matrix(points).transpose())
    B = q[:,:n]

    # compare numpy version with integer-value method--just pulling out a set of Linearly independent points
    ## get first nonzero point
    #point0 = points[0]
    #if np.count_nonzero(point0) < 1:
    #    for p in points:
    #        if np.count_nonzero(p) > 0:
    #            point0 = p
    #            break
    ## then find a LI set of points
    #LI_vecs = [point0]
    #current_rank = 1
    #for p in points:
    #    if np.count_nonzero(p) > 0:
    #        current_mat = np.matrix(LI_vecs + [p])
    #        if np.linalg.matrix_rank(current_mat) > current_rank:
    #            LI_vecs.append(p)
    #            current_rank = np.linalg.matrix_rank(current_mat)
    #B = np.matrix(LI_vecs).transpose()
    return B


def identicalLists(list1, list2):
    if len(list1) != len(list2):
        return False
    else:
        return np.all([list1[x] == list2[x] for x in range(len(list1))])


def intersectionDim(a, b):
    return np.linalg.matrix_rank(np.matrix((np.array(a)[:, None] == np.array(b)).all(-1).any(1)))


def intSolve(A, b, nonneg=False):
    sgn = "nonnegative" if nonneg else "none"
    return constrainSolve(A,b,sgn)

def constrainSolve(A, b, sgn="nonnegative"):
    nc = A.shape[1]
    x = cp.Variable(nc, integer=True)
    # set up the L2-norm minimization problem
    if sgn == "nonnegative":
        prob = cp.Problem(cp.Minimize(cp.norm(A @ x - b, 2)), [x >= 0])
    elif sgn == "nonpositive":
        prob = cp.Problem(cp.Minimize(cp.norm(A @ x - b, 2)), [x <= 0])
    else:
        prob = cp.Problem(cp.Minimize(cp.norm(A @ x - b, 2)))
    sol = prob.solve(solver = 'ECOS_BB')
    return np.array(x.value, dtype='int32').ravel()



def isAcyclic(mat):
    edges = edgesFromMatrix(mat)
    verts = range(mat.shape[0])

    for first_v in verts:
        visited = [0 if x != first_v else 1 for x in verts]
        result = findCycleDFS(first_v, visited, edges)
        if result:
            return False
    return True


def isClosedUnderArrows(V, Q):
    return np.all(Q[V,:].sum(axis=0) >= 0)


# checks if a graph(represented by node_list and edge_list)is connected
def isConnected(node_list, edge_list):
    disconnected = False
    p = node_list[0]
    for q in node_list[1:]:
        if not existsPathTo(p, q, edge_list):
            disconnected = True
            break
    return not disconnected



def matDiff(a, b):
    if a.shape != b.shape:
        return True
    return np.all(a != b)



def matrixFromEdges(edges, oriented=True):
    nv = len(set(np.array(edges).flatten()))
    tail = -1 if oriented else 1
    return np.matrix([[tail if x == e[0] else 1 if x == e[1] else 0 for x in range(nv)] for e in edges]).transpose()


def matSolve(a=None, x=None, b=None):
    # solves the system ax = b where all 3 are 
    # (possibly non-square) matrices, and the unknown 
    # to solve for is the argument that is not defined. 
    # e.g.: to solve the system AX = B for B, 
    # try matSolve(a=A, x=X) and similarly for A and X

    if (x is not None) and (b is not None):
        # solve for a by multiplying b by a right inverse of x
        xinv = np.linalg.pinv(x)
        return b*xinv
    elif (a is not None) and (b is not None):
        # solve for x by inverting a
        ainv = np.linalg.pinv(a)
        return ainv*b
    elif (a is not None) and (x is not None):
        return a*x


def primalUndirectedCycle(G):
    """gives the edges that comprise an undirected cycle in the graph G, 
       (which is assumed to contain a single cycle) and returns the ordered cycle

        input: G(list of tuples): edges of graph G
        output: cycle(list of tuples): tuple representation of the edges contained in the cycle
    """
    for i, edge in enumerate(G):
        is_cycle, cycle = existsPathTo(edge[1], edge[0], G[:i] + G[i + 1:], False, True, [edge])
        if is_cycle:

            # go through and index each edge so as to allow for multiple edges between vertices
            edge_indices = []
            met_edges = []
            for c_e in cycle:
                for g_i, g_e in enumerate(G):
                    if g_i not in met_edges :
                        if g_e == c_e: 
                            met_edges.append(g_i)
                            edge_indices.append(g_i)
                            break
                        elif (g_e[1], g_e[0]) == c_e:
                            met_edges.append(g_i)
                            edge_indices.append(-(g_i+1))
                            break
            return edge_indices
    return []


def proj(a, b):
    # projects vector a onto basis (consisting of a list b of vectors)
    return np.array(np.matrix([bb*np.dot(a, bb)/np.dot(bb, bb) for bb in b if np.dot(bb,bb) > 0]).sum(axis=0).tolist()[0])


def spanningTree(M, tree_format="edge"):
    """ returns the first spanning tree for the graph Q that is found
    """
    Q0,Q1 = M.shape

    all_edges = edgesFromMatrix(M)
    all_nodes = range(Q0)

    edge_indices = []
    
    d = Q1 - Q0 + 1
    if d > 0:
        # try removing every combination of d edges and see if result is a tree
        d_tuples_to_remove = combinations(range(Q1), d)
        edges_kept = []
        edges_removed = []

        for d_tuple in d_tuples_to_remove:

            edges_kept = [e for i, e in enumerate(all_edges) if i not in d_tuple]
            edges_removed = [all_edges[d] for d in d_tuple]

            if isConnected(all_nodes, edges_kept) and isAcyclic(matrixFromEdges(edges_kept)):
                if tree_format != "edge":
                    edges_kept = [i for i in range(Q1) if i not in d_tuple]
                    edges_removed = list(d_tuple)
                return edges_kept, edges_removed

    if tree_format == "edge":
        return all_edges, []
    else:
        return range(Q1),[]


def subsetsClosedUnderArrows(mat):
    current_vertices = range(mat.shape[0])
    for i in current_vertices[1:]:
        for c in combinations(current_vertices, i):
            if isClosedUnderArrows(c, mat):
                yield c
    #return [c for i in current_vertices[1:] for c in combinations(current_vertices, i) if isClosedUnderArrows(c, mat)]
