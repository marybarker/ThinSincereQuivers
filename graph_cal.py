import numpy as np
from itertools import product, combinations, permutations
from itertools import combinations_with_replacement as all_combos
import collections


def is_a_column_permutation_of(A, B):
    """ checks if two matrices A and B are identical up to permutation of columns
        inputs:
            - A(numpy matrix)
            - B(numpy matrix)
        outputs:
            - bool of statement "A is a column permutation of B"
    """
    if A.shape != B.shape:
        return False
    a,b = A.T.tolist(), B.T.tolist()
    d = collections.defaultdict(int)
    for x in a:
        d[tuple(x)] += 1
    for x in b:
        d[tuple(x)] -= 1
    return not any(d.values())#itervalues())


# take a list of graphs (matrices, really) and return the unique elements. 
# i.e. a generalization of set() for lists of matrices
def unique_up_to_isomorphism(list_of_graphs):
    if len(list_of_graphs) > 1:
        to_remove = []
        for i, gi in enumerate(list_of_graphs):
            for j in range(i+1, len(list_of_graphs)):
                gj = list_of_graphs[j]
                if np.all((gi == gj)) or np.all((gi == -gj)) or is_a_column_permutation_of(gi, gj):
                    to_remove.append(j)

        to_remove = list(set(to_remove))
        return [x for i, x in enumerate(list_of_graphs) if not (i in to_remove)]
    return list_of_graphs

# generate all graphs (unique up to isomorphism) that have the 
# appropriate number of edges and vertices for a given value d
def generate_graphs(d):
    connectivity_matrices = []
    # first generate G0:
    G0s = [x for x in range(1, 2*(d - 1) + 1)]
    for G0 in G0s:
        # generate G1:
        G1 = G0 + d - 1

        list_of_aij_possibilities = [[0, 1, 2] for x in range(G0)]
        all_possible_columns = [x for x in product(*list_of_aij_possibilities) if sum(x) == 2]
        list_of_col_indices = [list(range(len(all_possible_columns))) for x in range(G1)]
        all_col_combinations = list(product(*list_of_col_indices))

        As = [np.matrix(np.column_stack([all_possible_columns[j] for j in i])) for i in all_col_combinations]
        As = [A for A in As if np.all(np.array((A.sum(axis=1)) >= 3))]

        # now make sure the graphs are unique up to graph isomorphism
        As = unique_up_to_isomorphism(As)
        connectivity_matrices.extend(As)
    return connectivity_matrices


# This set of functions deals with part d in the theorem.
# namely, it checks that each edge is contained in a cycle.
def edges_out_of(p, edge_list, oriented=False):
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

def exists_path_to(p, q, edge_list):
    """ checks if there exists a path from vertex p to vertex q that 
        can be obtained by stringing together edges from edge_list. 
    """
    is_path = False
    for edge in edges_out_of(p, edge_list):
        i, e = edge
        # extract endpoints of edge
        v = e[0]
        if p == e[0]:
            v = e[1]

        # check if there's an edge between p & q
        if q == v:
             is_path = True
             break
        else:
            new_edge_list = edge_list[:i] + edge_list[i + 1:]
            if exists_path_to(v, q, new_edge_list):
                is_path = True
                break
    return is_path


def in_a_cycle(e_index, edge_list):
    """
       inputs: 
           - e_index(int): index of edge to check
           - edge_list(list of lists): list of edges, each of which is of the form [v1, v2]
       output: 
           - boolean of the statement "edge e_index is contained in a cycle"

       NOTE: if every edge in a graph G is contained 
       in a cycle, then this algorithm will return a value of 
       true for each edge e. If only some are contained in a cycle, 
       then this algorithm might fail to detect a cycle. 

       TLDR: this function only works for the case I wrote it for. Don't reuse without checking.
    """
    success = False

    e = edge_list[e_index]
    E = edge_list[:e_index] + edge_list[e_index + 1:]

    # get other vertex associated to the edge
    p,q = e

    # check if there is a path from p to q that does not contain edge e
    if exists_path_to(q, p, E):
        success = True
    return success


def all_edges_in_a_cycle(edge_list):
    return np.all([in_a_cycle(e, edge_list) for e in range(len(edge_list))])


# FOR ORIENTED GRAPH CYCLE CHECKING
def exists_cycle(edge_list, vertex_list):
    visited = [0 for x in vertex_list]
    v = edge_list[0][0]
    visited[v] = 1
    return dfs_find_cycle(v, visited, edge_list)


# FOR ORIENTED GRAPH CYCLE CHECKING
def dfs_find_cycle(starting_vertex, visited, edge_list):
    for e_index, e in edges_out_of(starting_vertex, edge_list, oriented=True):
        if visited[e[1]]:
            return True
        visited[e[1]] = 1
        return dfs_find_cycle(e[1], visited, edge_list)
    return False


# checks if a graph(represented by node_list and edge_list)is connected
def is_connected(node_list, edge_list):
    disconnected = False
    p = node_list[0]
    for q in node_list[1:]:
        if not exists_path_to(p, q, edge_list):
            disconnected = True
            break
    return not disconnected


# returns list of tuples of the form (node1, node2) that corresponds to 
# an edge between node1 and node2. Loops are represented as (node1,node1)
def edges_of_graph(A, oriented=False):
    node_list = range(A.shape[0])
    edges = np.asarray(A.T)
    edge_list = [tuple(i for i in node_list if e[i] != 0) for e in edges]
    edge_list = [x if len(x) > 1 else (x[0], x[0]) for x in edge_list]

    if oriented:
        # check that ordering of nodes is done correctly 
        # based on the +-1 values in A not on the natural ordering 0-#vertices
        for i, e in enumerate(edge_list):
            if A[e[1],i] < 0:
                edge_list[i] = (e[1], e[0])
    return edge_list



def spanning_tree(Q):
    """ Returns a spanning tree of the quiver Q with |Q_1| - |Q_0| + 1 edges removed. 
        NOTE: if such a spanning tree is not possible, then it returns empty lists

        input: 
            - Q(numpy matrix): Matrix representation of quiver
        outputs:
            - edges_kept(list of tuples): list of the edges in the spanning tree
            - edges_removed(list of tuples)): list of the edges in the complement of the spanning tree
    """

    # edges of quiver Q represented as a list of tuples
    all_edges = edges_of_graph(Q, True)
    all_nodes = range(Q.shape[0])

    # number of edges to remove from spanning tree
    d = Q.shape[1] - Q.shape[0] + 1

    d_tuples_to_remove = list(combinations(range(len(all_edges)), d))
    for d_tuple in d_tuples_to_remove:
        edges_kept = [e for i, e in enumerate(all_edges) if i not in d_tuple]
        edges_removed = [all_edges[d] for d in d_tuple]
        if is_connected(all_nodes, edges_kept):
            return edges_kept, edges_removed
    return [], []


def get_directed_path_to(p, q, edge_list, edges_added):
    """ checks if there exists a path from vertex p to vertex q that 
        can be obtained by stringing together edges from edge_list, 
        and returns that (re-oriented) list of edges. 
    """
    is_path = False
    edges_to_add = list(edges_added)
    for edge in edges_out_of(p, edge_list):
        i, e = edge
        # extract endpoints of edge
        coef = -1
        v = e[0]
        if p == e[0]:
            v = e[1]
            coef=1

        # check if there's an edge between p & q
        if q == v:
             is_path = True
             edges_to_add.append(coef*(i+1))#(p, q))
             break
        else:
            new_edge_list = edge_list[:i] + edge_list[i + 1:]
            a, b = get_directed_path_to(v, q, new_edge_list, edges_to_add + [coef*(i+1)])#[(p, v)])
            if a:
                is_path = True
                edges_to_add = b
                break
    return is_path, edges_to_add


def primal_undirected_cycle(G):
    """gives the edges that comprise an undirected cycle in the graph G, 
       (which is assumed to contain a single cycle) and returns the ordered cycle

        input: G(list of tuples): edges of graph G
        output: cycle(list of tuples): tuple representation of the edges contained in the cycle
    """
    for i, edge in enumerate(G):
        is_cycle, cycle = get_directed_path_to(edge[1], edge[0], G[:i] + G[i + 1:], [i+1])
        if is_cycle:
            return cycle
    return []



def flow_polytope(Q, W=None):
    """ return the vertices of the flow polytope associated to the quiver Q
    input: Q(numpy matrix): #vertices x #arrows matrix corresponding to quiver Q
    output: vertices: list of vertices comprising the convex hull of the flow polytope
    """

    T, edges_removed = spanning_tree(Q)
    edges = T + edges_removed

    # make sure the weights are the canonical ones. 
    if not W:
        W = [1 for x in range(len(edges))]

    f = []
    for i, edge in enumerate(edges_removed):
        cycle = primal_undirected_cycle(T+[edge])
        fi = np.zeros(len(edges))
        for j in cycle:
            if j > 0:
                fi[j-1] = W[j-1]
            elif j < 0:
                fi[1 - j] = -W[1 - j]
        f.append(fi)
    vertices = [list(f[i][j] for i in range(len(edges_removed))) for j in range(len(edges))]
    return vertices


def Step1(n):
    # obtain matrices for all graphs with d=n
    Temporary = [(M, edges_of_graph(M)) for M in generate_graphs(n)]
    return [y for (y, x) in Temporary if all_edges_in_a_cycle(x) and is_connected(range(y.shape[0]), x)]



def Step2(mat):
    # remove all loops from graph with connectivity matrix mat
    orig_num_edges = mat.shape[1]
    loops_broken = []

    A = mat.copy()
    B = np.isin(A, [2])
    if B.sum() > 0:
        pairs = []
        for r, rval in enumerate(B):
            for c, cval in enumerate(B[r]):
                if B[r][c]:
                    pairs.append((r,c))

        for pair in pairs:
            A[pair[0],pair[1]] = 1
            col_to_add = A[:, pair[1]]

            A = np.append(A, col_to_add, axis=1)
            row_to_add = np.zeros(A.shape[1])
            row_to_add[pair[1]] = 1
            row_to_add[-1] = 1

            A = np.append(A, np.matrix(row_to_add), axis=0)
            loops_broken.append(pair[1])
        loops_broken.extend([x for x in range(orig_num_edges, A.shape[1])])
    return A, loops_broken



def Step3(mat, edges_to_split):
    # split the corresponding list of edges in matrix mat
    # NOTE: assuming all entries are 0 and 1 (no loops, which would have value 2)
    A = mat.copy()
    nr,nc = A.shape
    for edge in edges_to_split:
        column = [i for i, x in enumerate(mat.T[edge,:].flatten().tolist()[0]) if x == 1]
        row1, row2 = column[0], column[1] # THIS WILL BREAK if mat is not in correct form
        A = np.insert(A, nc, 0, axis=1) # add column of 0s to matrix
        A = np.insert(A, nr, 0, axis=0) # add row of 0s
        A[row2, edge] = 0
        A[row2, nc] = 1
        A[nr, edge] = 1
        A[nr, nc] = 1
    return A



def Step4(mat):
    # put signs on the arrows.

    #first step: find the vertices with valence 2 (i.e. row sum = 2)
    rowsums = np.array(mat.T.sum(axis=0).flatten().tolist()[0])
    val2_verts = np.where(rowsums == 2)[0]
    diag = np.zeros((mat.shape[0], mat.shape[0]))
    np.fill_diagonal(diag, [1 if not (x in val2_verts) else -1 for x in range(mat.shape[1])])
    mat = np.matmul(diag, mat)
    mats = [mat]

    # now for each column that does not add up to 0, generate all possible
    # combinations of +-1 for it
    A = mat.T.copy()
    columns_to_alter = np.where(np.array(A.sum(axis=1).flatten().tolist()[0]) != 0)[0]
    columns_to_save = [x for x in range(A.shape[0]) if x not in columns_to_alter]

    if len(columns_to_alter) > 0:
        mats=[]
        possible_columns = []
        for col in columns_to_alter:
            column = A[col,:].tolist()[0]
            idx = column.index(1)

            column[idx] = -1
            col1 = column
            col2 = [-x for x in column]
            possible_columns.append((col1, col2))

        base_mat = mat[:, columns_to_save]
        col_choices = list(all_combos([0,1], len(possible_columns)))

        mats = [np.column_stack((base_mat, np.column_stack(A))) for A in [[possible_columns[i][y] for i, y in enumerate(x)] for x in col_choices]]

    #save only the graphs that are cycle-free
    mats = [-x for x in mats if not exists_cycle(edges_of_graph(x, True), range(x.shape[0]))]
    return unique_up_to_isomorphism(mats)



def Step5(mats):
    # take a list of quivers/directed graphs and return a list unique up to isomorphism.

    to_remove = np.zeros(len(mats))
    for iM, M in enumerate(mats):
        M_shape = M.shape
        M_valences = sorted(list(M.sum(axis=1)))

        # go through remaining matrices to see if there is a match
        for iN in range(iM+1, len(mats)):
            if not to_remove[iN]:
                N = mats[iN]

                if np.all(N == M):
                    # check if they're equal as matrices
                    to_remove[iN] = 1

                elif N.shape == M_shape:
                    # check if shapes are the same (least expensive)
                    N_valences = sorted(list(N.sum(axis=1)))

                    # then check graph valency the same
                    if N_valences == M_valences:

                        # finally see if one is a permutation of the vertices of the other
                        possible_permutations = list(permutations(range(N.shape[0])))
                        for possible in possible_permutations:
                            N_p = N[possible,:]
                            if is_a_column_permutation_of(N_p, M):
                                to_remove[iN] = 1
                                break
    return [m for i, m in enumerate(mats) if not to_remove[i]]
