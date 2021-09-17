import ThinSincereQuivers as tsq
import graph_tools as gt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


#a = [[0,2],[0,3],[0,4],[1,2],[1,3],[1,4]]
#
#Q1 = tsq.ToricQuiver(a)
#print(Q1)
#
#Q2 = tsq.ToricQuiver(tsq.matrixFromEdges(a), flow=[-1,2,-3,4,-5,6])
#print(Q2)
#
#Q3 = tsq.chainQuiver([1,2,3])
#print(Q3)
#
#print(tsq.allSpanningTrees(Q1, tree_format="vertex"))
#
#print("unstable subquivers are: ", list(tsq.unstableSubquivers(tsq.bipartiteQuiver(2,3,[2,2,1,0,0,1]), output_format="list")))
#
#print(tsq.basisForFlowPolytope(Q1,[0,1,2,3])) 
#
#print(tsq.isClosedUnderArrows([1,2,3], Q1))
#print(tsq.isClosedUnderArrows([2,3], Q1))
#
#print(tsq.flowPolytope(Q1, polytope_format="a"))
#print("\n")
#print("\n")
#print(tsq.flowPolytope(Q1))
#print("\n")
#
#print([x for x in tsq.subsetsClosedUnderArrows(Q1.incidence_matrix)])
#
#print(tsq.maximalUnstableSubquivers(Q1, True))
#
#print(tsq.maxCodimensionUnstable(Q1))
#print(list(tsq.primitiveArrows(Q1)))
#
#print(list(tsq.primitiveArrows(tsq.ToricQuiver([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]))))
#
#
#print(tsq.potentialWalls(Q1))
#
#print(tsq.makeTight(Q1, [-5,-1,2,2,2]))
#
#b = tsq.bipartiteQuiver(2,3,flow=[2,2,1,0,0,1])
#print(tsq.maximalUnstableSubquivers(b, True))
#print(list(tsq.stableTrees(b, b.weight)))

#Q3 = tsq.ToricQuiver([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
#cs = tsq.coneSystem(Q3)
#print(list(cs))
#print(len(list(cs)))


#print(tsq.mergeOnArrow(Q1, 5, Q1, 0))
#print(tsq.mergeOnVertex(Q1, 1, Q1, 0))

#print(len(list(tsq.coneSystem(Q1))))
#K5 = tsq.ToricQuiver([[0,1],[0,2],[0,3],[0,4],[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]])
#print(len(list(tsq.coneSystem(K5))))

def view_conesystem_3d():
    Q = tsq.ToricQuiver([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])

    spanning_trees = tsq.allSpanningTrees(Q, tree_format="vertex")
    qcmt = Q.incidence_matrix.transpose()

    Q_dim = qcmt.shape[0] - qcmt.shape[1] + 1
    cone_dim = qcmt.shape[1] - 1

    # create list of cones CT for each tree T
    tree_chambers = list(map(lambda st: qcmt[st[0],:], spanning_trees))
    primitive_arrows = tsq.primitiveArrows(Q)

    # vector pointing from centroid of cone to the origin
    vec_to_conebase = -(1./len(primitive_arrows))*np.matrix(qcmt[primitive_arrows,:]).sum(axis=0).transpose()
    # find Vertices in C(Q) that define the subchambers
    V = set([tuple(r) for mat in tree_chambers for r in mat.tolist()])
    V = [np.array(v) for v in V]

    # find a basis for the lower-dimensional space that is spanned by the vertices in C(Q)
    B = gt.findLowerDimSpace(V)
    A = np.dot(np.linalg.inv(np.dot(B.transpose(),B)),B.transpose())

    # convert all vertices to this lower space to define cones
    lower_dim_centroid = np.array(np.dot(A, -vec_to_conebase).transpose().tolist()[0])
    lower_dim_conebase = [-lower_dim_centroid]
    lower_dim_trees = [tsq.Cone( \
                              np.dot(t, A.transpose()).tolist()+lower_dim_conebase) \
                          for t in tree_chambers]

    # find all pairs of cones with full-dimensional interesection
    aij = [[i+j+1 \
        for j, tcj in enumerate(lower_dim_trees[i+1:]) \
        if (tci.intersection(tcj).dim >= cone_dim)] \
        for i, tci in enumerate(lower_dim_trees) \
    ]

    minv=100
    maxv=-100
    for i, tci in enumerate(lower_dim_trees):
        for j, tcj in enumerate(lower_dim_trees[i+1:]):
            print("trying to intersect: %d and %d"%(i, i+j+1))
            ij = tci.intersection(tcj, printout=True, B=B)
            print(tree_chambers[i])
            print(tree_chambers[i+j+1])
            minv=min(minv, np.amin(np.array(tci.vertices)))
            maxv=max(maxv, np.amax(np.array(tci.vertices)))

    # plot the lower dimensional version of the chamber
    ax = plt.axes(projection='3d')
    cb = np.array(lower_dim_conebase[0])
    cbdir = cb / (cb[0]**2 + cb[1]**2 + cb[2]**2)**0.5
    azim = np.round(180*np.arctan2(cbdir[1], cbdir[0])/np.pi)
    elev = np.round(180*np.arctan2(1, cbdir[2])/np.pi)

    ax.view_init(elev, azim)
    ax.scatter3D([cb[0]], [cb[1]], [cb[2]]) # point showign the base of the cone
    for i, tci in enumerate(lower_dim_trees):
        pts = np.array(tci.vertices).transpose().tolist()
        ax.plot_trisurf(pts[0], pts[1], pts[2])

        pts1 = [[pts[0][x] for x in [0,1,2,0,2,3,0,1,3]], \
                [pts[1][x] for x in [0,1,2,0,2,3,0,1,3]], \
                [pts[2][x] for x in [0,1,2,0,2,3,0,1,3]]]

        ax.plot3D(pts1[0], pts1[1], pts1[2])
        ax.set_xlim3d(minv,maxv)
        ax.set_ylim3d(minv,maxv)
        ax.set_zlim3d(minv,maxv)

        plt.draw()
        plt.pause(2)
    plt.pause(20)

view_conesystem_3d()
