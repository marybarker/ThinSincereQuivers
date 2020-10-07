needsPackage("ThinSincereQuivers")

Q = bipartiteQuiver(2, 3)
assert isTight Q
assert (theta Q === {-3,-3,2,2,2})
assert (sort entries transpose flowPolytope Q === {{-1, -1}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}, {1, 1}} )

th = {-5,-1,2,2,2}
F = incInverse(Q, th)
assert not isTight(Q, F)
assert (makeTight(Q, th) == toricQuiver({{0,1},{0,1},{0,1}}, {-1,1,1}))

Q = bipartiteQuiver(2, 2)
assert (subquivers(Q, Format=>"list") === {{0}, {1}, {2}, {3}, {0, 1}, {0, 2}, {1, 2}, {0, 3}, {1, 3}, {2, 3}, {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}})
