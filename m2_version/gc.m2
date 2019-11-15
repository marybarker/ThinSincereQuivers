loadPackage "Graphs"
------------------------------------------------------------
aslist = x -> (
    if class(x) === List then(
        x
    )
    else if class(x) === Sequence then(
        toList(x)
    )
    else if class(x) === Set then(
        toList(x)
    )
    else
        {x}
)
------------------------------------------------------------

------------------------------------------------------------
-- add all elements of a list together -- 
sumlist = {Axis => "None"} >> opts -> x -> (
    if opts.Axis == "None" then (
        s = sum(aslist(x));
    )
    else if opts.Axis == "row" then (
        s = flatten(for i in x list(sumlist(i)));
    )
    else if opts.Axis == "col" then (
       pivoted = entries(transpose(matrix(x)));
       s = flatten(for i in pivoted list(sumlist(i)));
    );
    s
)
------------------------------------------------------------


------------------------------------------------------------
-- take all possible combinations of length k from list l -- 
-- optional arguments: 
-- -- R(true/false) = with replacement
-- -- Minsum(numeric value) = exclude all combinations with sum below Minsum
-- -- Maxsum(numeric value) = exclude all combinations with sum above Maxsum
combinations = {R => true, Minsum => -1000, Maxsum => -1000} >> opts -> (l, k) -> (
    if k > 1 then (
        -- if we are using combinations with replacement -- 
        if opts.R then (
           mylist = flatten(join(for i in l list(for j from 0 to k - 1 list(i))));
           combs1 = unique(subsets(mylist, k));
           combs2 = unique(subsets(mylist, k));
           for i in combs2 do (combs1 = append(combs1, reverse(i)));
           combs = unique(combs1);
        )
        else (
           mylist = flatten(for i in l list(i));
           combs1 = unique(subsets(mylist, k));
           combs2 = unique(subsets(mylist, k));
           for i in combs2 do (combs1 = append(combs1, reverse(i)));
           combs = unique(combs1);
        );
    )
    else
        combs = for i in l list(i);

    -- if we are using restricting either by a minimum or maximum sum -- 
    if opts.Minsum != -1000 then (
       combs = for i in combs list(if sumlist(i) < opts.Minsum then (continue;) else (i));
    );
    if opts.Maxsum != -1000 then (
       combs = for i in combs list(if sumlist(i) > opts.Maxsum then (continue;) else (i));
    );
    combs
)
------------------------------------------------------------


------------------------------------------------------------
isPerm = (x, y) -> (
    toRet = "False";
    ax = entries(transpose(x));
    ay = entries(transpose(y));
    vals = toList (0..#ax - 1);
    allPermutations = permutations(vals);
    for perm in allPermutations do (
        if ax_perm == ay then (
            toRet = "True";
            break;
        );
    );
    toRet
)
------------------------------------------------------------


------------------------------------------------------------
-- get unique entries from list of undirected graphs -- 
unorientedUniqueUpToPermutation = x -> (
    if #x > 1 then (
        toSave = (0..(#x - 1));
        for i from 0 to #x - 2 do (
            for j from i + 1 to #x - 1 do (
                if isPerm(x#i, x#j) == "False" then (
                    continue;
                )
                else (
                    toSave = delete(j, toSave);
                )
            )
        );
        for i in toSave list(x#i)
    )
    else x
)
------------------------------------------------------------


------------------------------------------------------------
-- create all possible undirected graphs with #vertices=x#0 and #edges=x#1 -- 
allPossibleBaseGraphsForPair = (x) -> (
   g0 = x#0;
   g1 = x#1;

   -- get all possible columns for connectivity matrix (entries must be 0,1, or 2)
   possibleCols = combinations({0,1,2}, g0, R=>true, Minsum=>2, Maxsum=>2);

   -- all combinations of columns to create rows
   rowCombs = combinations((0..(#possibleCols - 1)), g1, R=>true);
   candidateMats = for i in rowCombs list(for j in i list(aslist(possibleCols#j)));
   candidateMats = for i in candidateMats list(transpose(matrix(i)));
   candidateMats = for i in candidateMats list(if min(sumlist(entries(i), Axis=>"row")) >= 3 then (i) else continue);
   candidateMats = unorientedUniqueUpToPermutation(candidateMats);
   for m in candidateMats list(if graphIsConnected(m) then (m) else (continue;))
)
------------------------------------------------------------


------------------------------------------------------------
-- create all possible undirected base graphs for quiver generation in dimension d-- 
undirectedGraphs = (d) -> (
   gPairs = apply((1..2*(d - 1)), i -> {i, i+d-1});
   connectivityMatrices = allPossibleBaseGraphsForPair \ gPairs;
   connectivityMatrices
)
------------------------------------------------------------


------------------------------------------------------------
-- yield the edges of a graph in the form of a list of pairs 
-- (v1, v2), where edge E is from v1 to v2
graphEdges = {Oriented => "False", RavelLoops => false} >> opts -> (G) -> (
    if opts.Oriented == "True" then (
        E = for e in entries(transpose(G)) list(
            {position(e, i -> i < 0), 
             position(e, i -> i > 0)}
        );
    )
    else (
        E = for e in entries(transpose(G)) list(
            positions(e, i -> (i > 0 or i < 0))
        );
        if opts.RavelLoops == true then (
            E = for e in E list(if #e > 1 then (e) else ({2:e#0}));
        );
    );
    E
)
------------------------------------------------------------


------------------------------------------------------------
-- check if graph is connected
graphIsConnected = (G) -> (
    gEdges = graphEdges(G);
    if max(for e in gEdges list(#e)) > 1 then (
        isConnected(graph(gEdges))
    )
    else (
        if max(flatten(gEdges)) < 1 then (
            true
        )
        else (
            false
        )
    )
)
------------------------------------------------------------


edgesOutOfPoint = {Oriented => false} >> opts -> (p, E) -> (
    if opts.Oriented then (
        for i from 0 to #E - 1 list(e = E#i; if p != e#0 then (continue;) else (i, e))
    )
    else (
        for i from 0 to #E - 1 list(e = E#i; if (p != e#0 and p != e#1) then (continue;) else (i, e))
    );
)
------------------------------------------------------------
-- check if there exists a path between p and q by appending 
-- edges in E(which is a list of pairs (v1, v2). 
-- optional arguments: 
-- -- Oriented(true/false) = whether or not the graph should be oriented
-- -- savePath(true/false) = whether or not to return the edges involved in the path
-- -- edgesAdded(list) = internal mechanism for computing for savePath
pathBetween = {Oriented => false, savePath => false, edgesAdded => {}} >> opts -> (p, q, E) -> (

    existsPath = false;
    currentPath = {};
    for edge in edgesOutOfPoint(p, E, Oriented => opts.Oriented) do (

        --- get the edge index and enpoints
        i = edge#0;
        e = edge#1;
        v = e#1;
        if p == e#1 then (
            v = e#0;
        );
        if opts.savePath then (
            currentEdges = append(list(opts.edgesAdded), {p, v})
        );

        if q == v then (
            existsPath = true;
            break;
        )
        else (
            remainingEdges = for j from 0 to #E - 1 list(if j == i then (continue;) else E#j);
            if pathBetween(v, q remainingEdges, Oriented => opts.Oriented, savePath => opts.savePath, edgesAdded=>currentEdges) then (
                existsPath = true;
                break;
            )
        );
    );
    if opts.savePath then (
        (existsPath, currentPath)
    )
    else (
        existsPath
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- optionally if NautyGraphs generateGraphs() were working, we could do: 
---- 
---- undirectedGraphs = (d) -> (
----    gPairs = apply((1..2*(d - 1)), i -> {i, i+d-1});
----    connectivityMatrices = apply(gPairs, i -> generateGraphs(n, e, MinDegree=>3, OnlyConnected=>true));
---- )
-- but unfortunately, it is not. 
------------------------------------------------------------
--- to run this in the M2 terminal, type: 
-- > input "gc.m2"
--- Now you have all of the functions in this module!!!
--- so, for example to get all possible graphs with 2 vertices and 4 edges, type: 
-- > gs = allPossibleBaseGraphsForPair({2,4})
