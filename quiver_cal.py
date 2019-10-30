import numpy as np
from itertools import combinations


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
    return [-x for x in M.sum(axis=1).flatten().tolist()[0]]



def is_stable(M, subM): # CHECK!! (with paper)
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
    
    return not np.any(sums <= 0)
    "if sum>0 for some closed subset. Then np.any(sums >0) = true"



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
    return [ subQ for subQ in L if not is_stable(M, list(subQ))]



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








QT=np.array([
  [ 1, 1, 1, 0, 0, 0],
  [-1, 0, 0,-1, 0, 0],
  [ 0,-1, 0, 0,-1, 0],
  [ 0, 0,-1, 0, 0,-1],
  [ 0, 0, 0, 1, 1, 1]
])

WQT=np.array([
  [ 2, 3, 4, 0, 0],
  [-2, 0, 0,-5, 0],
  [ 0,-3, 0, 0,-6],
  [ 0, 0,-4, 0, 0],
  [ 0, 0, 0, 5, 6]
])

Q = np.asmatrix(QT)
WQ = np.asmatrix(WQT)

len(subquivers(Q))
theta(WQ)
subsets_closed_under_arrows(WQ)
subsets_closed_under_arrows(Q)
a = is_stable(Q, [0,3,4,5])
b = is_stable(Q, [1,2,4,5])
print(a, b)


unstable_subquivers = all_unstable(Q)
len(unstable_subquivers)
A = all_maximal_unstable(Q)

for e in A:
    print(e)
    print("")


WT=np.array([
  [ 100, 100, 100,   0,   0,   0],
  [-100,   0,   0,  -1,   0,   0],
  [   0,-100,   0,   0,  -1,   0],
  [   0,   0,-100,   0,   0,  -1],
  [   0,   0,   0,   1,   1,   1]
])
W= np.asmatrix(WT)

is_tight(Q)

not_tight = -np.matrix([
     [ 1, 1, 0, 0, 0],
     [ 0, 0,-1, 0, 0],
     [-1, 0, 0,-1, 0],
     [ 0,-1, 0, 0,-1],
     [ 0, 0, 1, 1, 1]
])
all_maximal_unstable(not_tight)
is_tight(not_tight)
