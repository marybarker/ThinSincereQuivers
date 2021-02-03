needsPackage("ThinSincereQuivers")

Q = bipartiteQuiver(2, 3)
assert isTight Q
assert (theta Q === {-3,-3,2,2,2})
assert (flowPolytope Q == {{0, 0}, {0, -1}, {2, -2}, {1, -2}, {2, -1}, {1, 0}})
assert (entries basisForFlowPolytope Q === {{-1, 0}, {0, -1}, {1, 1}, {1, 0}, {0, 1}, {-1, -1}})

th = {-5,-1,2,2,2}
F = incInverse(Q, th)
assert not isTight(Q, F)
assert isAcyclic makeTight(Q, th)
assert (makeTight(Q, th) == toricQuiver({{1,0},{1,0},{1,0}}, {-1,1,1}))
assert not isSemistable(Q, {1,2})
assert isStable(Q, {1,2,3,4})
assert (stableTrees(th,Q) === {{0, 1, 2, 5}, {0, 1, 2, 4}, {0, 1, 2, 3}})


Q = bipartiteQuiver(2, 2)
assert (subquivers(Q, Format=>"list") === {{0}, {1}, {2}, {3}, {0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, {2, 3}, {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}})
assert (allSpanningTrees Q == {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}})
assert not isClosedUnderArrows({1, 2}, Q)

munsbs = maximalUnstableSubquivers Q
vals = flatten for key in keys(munsbs) list(munsbs#key)
assert (vals ===  {{0}, {1}, {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}})

P = chainQuiver {2,2}
Q = chainQuiver {1,2,3}

assert (mergeOnVertex(P,2,Q,0) == chainQuiver {2,2,1,2,3})
assert (mergeOnArrow(P,3,Q,0) == chainQuiver {2,2,2,3})
assert (neighborliness P === 2)
