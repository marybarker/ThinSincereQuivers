import numpy as np
import collections
from itertools import product, combinations, permutations
from itertools import combinations_with_replacement as all_combos


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
    return not any(d.values())


# take a list of graphs (matrices, really) and return the unique elements. 
# i.e. a generalization of set() for lists of matrices
def unoriented_unique_up_to_isomorphism(list_of_graphs):
    if len(list_of_graphs) > 1:
        to_remove = np.zeros(len(list_of_graphs))
        for i, gi in enumerate(list_of_graphs):
            for j in range(i+1, len(list_of_graphs)):
                if to_remove[j]:
                    pass
                else:
                    gj = list_of_graphs[j]
                    if np.all((gi == gj)) or is_a_column_permutation_of(gi, gj):
                        to_remove[j] = 1

        return [x for i, x in enumerate(list_of_graphs) if not (to_remove[i])]
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
        all_col_combinations = all_combos(range(len(all_possible_columns)), G1)
        As = [np.matrix(np.column_stack([all_possible_columns[j] for j in i])) for i in all_col_combinations]
        As = [A for A in As if np.all(np.array((A.sum(axis=1)) >= 3))]

        # now make sure the graphs are unique up to graph isomorphism
        As = unoriented_unique_up_to_isomorphism(As)
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


def exists_path_to(p, q, edge_list, oriented=False, returnPath=False, savedPath=[]):
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

    for edge in edges_out_of(p, edge_list, oriented=oriented):
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
                success, path = exists_path_to(v, q, remaining_edges, oriented, returnPath, savedPath)
            else:
                success = exists_path_to(v, q, remaining_edges, oriented, False)
            if success:
                is_path = True
                currentPath += path
                break
    if returnPath:
        return is_path, currentPath
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
def graph_from_edges(Es):
    all_nodes = list(set([y for x in Es for y in x]))
    num_edges = len(Es)
    A = np.zeros((len(all_nodes), num_edges))

    for i, e in enumerate(Es):
        A[e[0], i] = -1
        A[e[1], i] = 1
    return A


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
    """ Returns a spanning tree(the first one that is encountered) of 
        the quiver Q with |Q_1| - |Q_0| + 1 edges removed. 
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


def primal_undirected_cycle(G):
    """gives the edges that comprise an undirected cycle in the graph G, 
       (which is assumed to contain a single cycle) and returns the ordered cycle

        input: G(list of tuples): edges of graph G
        output: cycle(list of tuples): tuple representation of the edges contained in the cycle
    """
    for i, edge in enumerate(G):
        is_cycle, cycle = exists_path_to(edge[1], edge[0], G[:i] + G[i + 1:], False, True, [edge])
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


def flow_polytope(Q, W=None):
    """ return the vertices of the dual of the flow polytope associated to the quiver Q
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
        cycle = primal_undirected_cycle(T + [edge])

        # now replace the indices in cycle that are for edge so as to match indexing on list(edges)
        #to_replace = {len(T):len(T)+i, -(len(T)+1):-(len(T)+i+1)}
        cycle = [x + i if x == len(T) else -(len(T) + i + 1) if x == -(len(T) + 1) else x for x in cycle]

        fi = np.zeros(len(edges))
        for j in cycle:
            if j >= 0:
                fi[j] = W[j]
            elif j < 0:
                fi[-(1 + j)] = -W[-(1 + j)]
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
            A[pair[0], pair[1]] = 1
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
    return unoriented_unique_up_to_isomorphism(mats)



def Step5(mats):
    # take a list of quivers/directed graphs and return a list unique up to isomorphism.

    to_remove = np.zeros(len(mats))
    for iM, M in enumerate(mats):
        M_shape = M.shape
        M_valences = sorted(list(abs(M).sum(axis=1)))

        # go through remaining matrices to see if there is a match
        for iN in range(iM+1, len(mats)):
            if not to_remove[iN]:
                N = mats[iN]

                if np.all(N == M):
                    # check if they're equal as matrices
                    to_remove[iN] = 1

                elif N.shape == M_shape:
                    # check if shapes are the same (least expensive)
                    N_valences = sorted(list(abs(N).sum(axis=1)))

                    # then check graph valency the same
                    if N_valences == M_valences:

                        # finally see if one is a permutation of the vertices of the other
                        for possible in permutations(range(N.shape[0])):
                            N_p = N[possible,:]
                            if is_a_column_permutation_of(N_p, M):
                                to_remove[iN] = 1
                                break
    return [m for i, m in enumerate(mats) if not to_remove[i]]


def acyclic_quivers(d):
    graphs = Step1(d)

    step2_graphs, loops_broken = zip(*[Step2(A) for A in graphs])

    step3_graphs = []
    for i, A in enumerate(step2_graphs):
        n_edges = A.shape[1]
        edges_to_break = [x for x in range(n_edges) if not (x in loops_broken[i])]
        if len(edges_to_break) > 0:
            for y in range(len(edges_to_break)+1):
                for edge_list in combinations(edges_to_break, y):
                    step3_graphs.append(Step3(A, list(edge_list)))
        else:
            step3_graphs.append(A)
    step3_graphs = unoriented_unique_up_to_isomorphism(step3_graphs)

    step4_graphs = []
    for A in step3_graphs:
        step4_graphs.extend(Step4(A))
    step4_graphs = unoriented_unique_up_to_isomorphism(step4_graphs)

    step5_graphs = Step5(step4_graphs)

    return step5_graphs


def subquivers(M):
    """
    input: M(numpy matrix): connectivity matrix of some quiver Q
    output: L(list of matrices): all subquivers of M with same vertices
    """
    num_arrows = M.shape[1]
    arrows = range(num_arrows)
    L = [M[:,j] for i in arrows[1:] for j in combinations(arrows, i)]
    return L



def subsets_closed_under_arrows(M):
    """
    v.1: sum matrix columns

    input: M(numpy matrix): connetivity matrix of some quiver Q

    output: list of tuples(corresponding to vertices)
            all subsets of the vertices of Q that are closed under arrows.
            i.e. S subset Q0 such that for v in S, the arrows in Q with
            tail v have head in S
    """
    num_vertices = M.shape[0]
    vertices = range(num_vertices)
    return [j for i in range(1, num_vertices) for j in combinations(vertices, i) if np.all(M[j,:].sum(axis=0) <= 0)]



def theta(M):
    """
    input: M(numpy matrix): weighted connectivity matrix of some quiver Q

    output: theta(list): ordered list of weights for the vertices of Q
    """
    return [x for x in M.sum(axis=1).flatten().tolist()[0]]



def is_stable(M, subM):
    """
    inputs: M(numpy matrix): weighted connectivity matrix of some quiver Q
            subM(list): subset of the indices[0,...,#columns of M] corresponding to the arrows of the subquiver

    output: boolean of statement "M is W-stable"
    """
    vertices = list(set(r for r in range(M.shape[0]) if np.any(np.array([M[r,c] != 0 for c in subM])))) #of the subquiver
    theta_calc = theta(M)
    weights = [theta_calc[i] for i in vertices]

    submatrix = M[np.array(vertices),:]
    submatrix = submatrix[:,subM]
    # calculate the weight sums for each subset of M that is closed under arrows
    sums = np.array([sum(np.array(weights)[np.array(S)]) for S in subsets_closed_under_arrows(submatrix)])
    return np.all(sums < 0)



def all_unstable(M):
    """
    inputs: M(numpy matrix): connectivity matrix of some quiver Q
            W(list or numpy array): list of weights

    output: list of all unstable (wrt W) subquivers of M that are closed under arrows
    NOTE: unstable in this context includes things that are not stable i.e. unstable and semi-stable.
    """
    num_arrows = M.shape[1]
    arrows = range(num_arrows)
    L = [list(j) for i in arrows[1:] for j in combinations(arrows, i)]
    return [subQ for subQ in L if not is_stable(M, list(subQ))]



def contained_in(Q1, Q2):
    """check if Q1 is a subset of Q2
        inputs: two lists of vertices(interpreted as subquivers)
        output: boolean evaluating the statement Q1 subset of Q2
    """
    if np.array_equal(Q1, Q2):
        return False
    if len([x for x in Q1 if x in Q2]) == len(Q1):
        return True
    return False


def all_maximal_unstable(M):
    """
    """
    # all unstable subquivers of M closed under arrows
    # NOTE: unstable in this context includes things that are not stable i.e. unstable and semi-stable.
    unstable_list = all_unstable(M)

    answer = []
    for Q in unstable_list:
        maximal = True
        for Q_other in unstable_list:
            if contained_in(Q, Q_other):
                maximal = False
        if maximal:
            answer.append(Q)
    return answer


def is_tight(Q):
    """
    input: Q(weighted connectivity matrix for quiver)
    output: boolean of statement "Q is tight"
    """
    n_arrows = Q.shape[1]
    return np.all([len(x) != (n_arrows - 1) for x in all_maximal_unstable(Q)])


def neighborliness(Q):
    max_unstables = all_maximal_unstable(Q)
    n = [len(m) for m in max_unstables] + [0]
    return Q.shape[1] - max(n)
