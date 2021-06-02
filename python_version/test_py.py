import ThinSincereQuivers as tsq


a = [[0,2],[0,3],[0,4],[1,2],[1,3],[1,4]]

Q1 = tsq.ToricQuiver(a)
print(Q1)

Q2 = tsq.ToricQuiver(tsq.matrixFromEdges(a), flow=[-1,2,-3,4,-5,6])
print(Q2)

Q3 = tsq.chainQuiver([1,2,3])
print(Q3)

print(tsq.allSpanningTrees(Q1, tree_format="vertex"))

print("unstable subquivers are: ", list(tsq.unstableSubquivers(tsq.bipartiteQuiver(2,3,[2,2,1,0,0,1]), output_format="list")))

print(tsq.basisForFlowPolytope(Q1,[0,1,2,3])) 

print(tsq.isClosedUnderArrows([1,2,3], Q1))
print(tsq.isClosedUnderArrows([2,3], Q1))

print(tsq.flowPolytope(Q1, polytope_format="a"))
print("\n")
print("\n")
print(tsq.flowPolytope(Q1))
print("\n")

print([x for x in tsq.subsetsClosedUnderArrows(Q1.connectivity_matrix)])

print(tsq.maximalUnstableSubquivers(Q1, True))

print(tsq.maxCodimensionUnstable(Q1))
print(list(tsq.primitiveArrows(Q1)))

print(list(tsq.primitiveArrows(tsq.ToricQuiver([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]))))


print(tsq.potentialWalls(Q1, tsq.theta(Q1)))

print(tsq.makeTight(Q1, [-5,-1,2,2,2]))

b = tsq.bipartiteQuiver(2,3,flow=[2,2,1,0,0,1])
print(tsq.maximalUnstableSubquivers(b, True))
print(list(tsq.stableTrees(b, b.weight)))

Q3 = tsq.ToricQuiver([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
print(list(tsq.coneSystem(Q3)))
print(len(list(tsq.coneSystem(Q3))))
