restart
loadPackage("ToricQuiver")

Q = bipartiteQuiver(2, 3)

isTight Q
theta Q

subquivers Q
subquivers(Q, Format=>"list")

sq = first subquivers(Q, Format=>"list")

Q_sq
Q^sq

flowPolytope Q

th = {-5,-1,2,2,2}
F = flatten entries first incInverse(Q, th)
isTight(Q, F)

makeTight(Q, th)
