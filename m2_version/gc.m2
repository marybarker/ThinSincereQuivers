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
-- -- Order(true/false) = whether or not the ordering of combination values matters
combinations = {R => true, Minsum => -1000, Maxsum => -1000, Order => true} >> opts -> (l, k) -> (
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
        combs = for i in l list(aslist(i));

    -- if we are using restricting either by a minimum or maximum sum -- 
    if opts.Minsum != -1000 then (
       combs = for i in combs list(if sumlist(i) < opts.Minsum then (continue;) else (i));
    );
    if opts.Maxsum != -1000 then (
       combs = for i in combs list(if sumlist(i) > opts.Maxsum then (continue;) else (i));
    );
    if opts.Order != true then (
        combs = unique(
            for i in combs list(sort(i))
        );
    );
    combs
)
------------------------------------------------------------


------------------------------------------------------------
-- return the indices of the list l in order the values occur 
-- in the sorted list sort(l)
sortedIndices = (l) -> (
    sortedVals = unique(sort(l));
    flatten(for i in sortedVals list(positions(l, x -> x == i)))
)
------------------------------------------------------------


------------------------------------------------------------
isPerm = (x, y) -> (
    toRet = false;

    xrows = entries(x);
    yrows = sort(entries(y));
    xcols = entries(transpose(x));
    ycols = entries(transpose(y));

    if #xrows == #yrows and #xcols == #ycols then (
        rs = toList(0..#xrows - 1);
        rowPermutations = permutations(rs);
        for rPerm in rowPermutations do (
            if xrows_rPerm == yrows then (
                toRet = true;
                break;
            )
            else (
                xrowsP = matrix(xrows_rPerm);
                cs = toList(0..#xcols - 1);
                colPermutations = permutations(cs);
                for cPerm in colPermutations do (
                    if sort(entries(xrowsP_cPerm)) == yrows then (
                        toRet = true;
                        break;
                    );
                );
            );
        );
    );
    toRet
)
------------------------------------------------------------
------------------------------------------------------------


------------------------------------------------------------
-- get unique entries from list of undirected graphs -- 
unorientedUniqueUpToPermutation = x -> (
    if #x > 1 then (
        toSave = (0..(#x - 1));
        for i from 0 to #x - 2 do (
            for j from i + 1 to #x - 1 do (
                if #positions(toSave, a -> a == j) > 0 then (
                    if isPerm(x#i, x#j) == false then (
                        continue;
                    )
                    else (
                        toSave = delete(j, toSave);
                    )
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
   flatten(connectivityMatrices)
)
------------------------------------------------------------


------------------------------------------------------------
-- yield the edges of a graph in the form of a list of pairs 
-- (v1, v2), where edge E is from v1 to v2
graphEdges = {Oriented => false, RavelLoops => false} >> opts -> (G) -> (
    if opts.Oriented == true then (
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
            E = for e in E list(if #e > 1 then e else toList(2:e#0));
        );
    );
    E
)
------------------------------------------------------------


replaceInList = (i, v, l) -> (
    insert(i, v, drop(l, {i,i}))
)

------------------------------------------------------------
-- yield the matrix rep of graph, given a list of edges as ordered 
-- pairs (this is the opposite of graphEdges() function. 
graphFromEdges = {Oriented => false} >> opts -> E -> (
    -- first, if oriented graph, then make sure this is reflected. 
    tailVal = 1;
    if opts.Oriented == true then (
        tailVal = -1;
    );

    nVerts = max(flatten(E))+1;
    cols = for i in E list(
        row = (nVerts:0);
        aslist(replaceInList(i#0, tailVal, replaceInList(i#1, 1, row)))
    );
    transpose(matrix(cols))
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


------------------------------------------------------------
edgesOutOfPoint = {Oriented => false} >> opts -> (p, E) -> (
    if opts.Oriented then (
        for i from 0 to #E - 1 list(e = E#i; if p != e#0 then (continue;) else (i, e))
    )
    else (
        for i from 0 to #E - 1 list(e = E#i; if (#positions(e, j -> j == p) < 1) then (continue;) else (i, e))
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- check if there exists a path between p and q by appending 
-- edges in E(which is a list of pairs (v1, v2). 
-- optional arguments: 
-- -- Oriented(true/false) = whether or not the graph should be oriented
-- -- SavePath(true/false) = whether or not to return the edges involved in the path
-- -- EdgesAdded(list) = internal mechanism for computing for SavePath
pathBetween = {Oriented => false, SavePath => false, EdgesAdded => {}} >> opts -> (p, q, E) -> (
    existsPath = false;
    currentPath = {};
    pathsToSee = edgesOutOfPoint(p, E, Oriented => opts.Oriented);

    for edge in pathsToSee do (
        --- get the edge index and enpoints
        i = edge#0;
        e = edge#1;
        v = e#1;
        if p == e#1 then (
            v = e#0;
        );
        if opts.SavePath then (
            currentEdges = append(toList(opts.EdgesAdded), {p, v});
        );

        if q == v then (
            existsPath = true;
            break;
        )
        else (
            remainingEdges = for j from 0 to #E - 1 list(if j == i then (continue;) else E#j);

            if opts.SavePath then (
                (ifPath, Path) = pathBetween(v, q, remainingEdges, Oriented => opts.Oriented, SavePath => true, EdgesAdded => currentEdges);
            )
            else (
                ifPath = pathBetween(v, q, remainingEdges, Oriented => opts.Oriented, EdgesAdded => currentEdges);
            );
            if ifPath then (
                existsPath = true;
                break;
            );
        );
    );
    if opts.SavePath then (
        (existsPath, currentPath)
    )
    else (
        existsPath
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- checks if there is a cycle containing given edge. 
edgeInCycle = (i, E) -> (
    if #E > 1 then (
        e = E#i;
        if #e > 1 then (
            p = e#0;
            q = e#1;
        )
        else (
            p = e#0;
            q = e#0;
        );
        indicesToSave = drop(toList(0..(#E-1)), {i,i});
        pathBetween(p, q, E_indicesToSave)
    )
    else (
        false
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- create all addmissable base graphs for a given dimension n
Step1 = (n) -> (
    for m in undirectedGraphs(n) list(
        es = graphEdges(m, RavelLoops => true);
        if #es > 0 then (
            inCycle = for i from 0 to #es - 1 list( edgeInCycle(i, es) );
            if all(inCycle, i -> i == true) then (m) else (continue;)
        )
    )
)
------------------------------------------------------------


------------------------------------------------------------
splitLoops = m -> (
    Es = graphEdges(m);
    nVerts = #entries(m);
    loopsBroken = for i from 0 to #Es - 1 list(
        e = Es#i;
        if (#e < 2) or #(delete(e#0, e)) < 1 then (
            i
        )
        else (continue;)
    );
    altLB = flatten(for i in loopsBroken list(i));
    newEdges = flatten(for i from 0 to #Es - 1 list (
        e = Es#i;
        if #loopsBroken != #delete(i, loopsBroken) then (
            p = position(loopsBroken, x -> x == i);
            altLB = append(replaceInList(p, loopsBroken_p + p, altLB), loopsBroken_p + p + 1);
            {{e#0, nVerts+p}, {nVerts+p, e#0}}
        )
        else (
            {e}
        )
    ));
    m = graphFromEdges(newEdges);
    m, altLB
)
------------------------------------------------------------


------------------------------------------------------------
splitEdges = (m, E) -> (
    Es = graphEdges(m);
    nVerts = #entries(m);

    for i from 0 to #E - 1 do (
        ei = E#i;
        e = Es#ei;
        Es = append(replace(i, {e#0, nVerts+i}, Es), {nVerts+i, e#1});
    );
    graphFromEdges(Es)
)
------------------------------------------------------------


------------------------------------------------------------
-- remove all loops from a given graph given in matrix form
-- and return the set of all loop-free copies with every 
-- combination of the remaining edges split
Step2 = (m) -> (
    origNumEdges = #graphEdges(m);
    (M, loops) = splitLoops(m);
    originalLooplessEdges = set(0..origNumEdges - 1) - set(loops);

    if #originalLooplessEdges > 0 then (
        otherEdgeCombs = for i from 1 to #originalLooplessEdges list(
            aslist(combinations(aslist(originalLooplessEdges), i, R => false, Order => false))
        );
        outputs = flatten(for i from 0 to #originalLooplessEdges - 1 list(
            for comb in otherEdgeCombs#i list(
                splitEdges(M, aslist(comb))
            )
        ));
        outputs = append(outputs, M);
        unorientedUniqueUpToPermutation(outputs)
    )
    else (
        {M}
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- put set of all admissable orientations on the matrix m
Step4 = (m) -> (
    rows = entries(m);

    -- first make all vertices of valence 2 into sinks. 
    m = matrix(for i in rows list(
        if sumlist(i) == 2 then (-1*i) else (i)
    ));
    cols = entries(transpose(m));
    es = graphEdges(m);

    -- now take every combination of +-1 for the columns that do 
    -- not yet add up to 0
    columnsToChange = for i from 0 to #cols - 1 list(
        if sumlist(cols#i) == 0 then (continue;) else (i)
    );

    if #columnsToChange > 0 then (
        columnChoices = for i in columnsToChange list(
            edge = es#i;
            tail = edge#0;
            head = edge#1;
            {replaceInList(tail, -1, cols#i), replaceInList(head, -1, cols#i)}
        );

        combinationsOfChoices = combinations({0, 1}, #columnsToChange, R=>true, Order=>false);
        for choice in combinationsOfChoices list(-1*transpose(matrix(
            for c from 0 to #cols - 1 list(
                if #columnsToChange == #delete(c, columnsToChange) then (
                    cols#c
                )
                else (
                    p = position(columnsToChange,i -> i == c);
                    colOpts = columnChoices#p;
                    choiceForCol = choice#p;
                    colOpts#choiceForCol
                )
            )
        )))
    )
    else (
        {-1*m}
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- take a list of quivers(directed graphs) and return the 
-- unique up to isomorphism list 
Step5 = (M) -> (
    if #M > 1 then (
        toSave = aslist(0..(#M - 1));
        for mi from 0 to #M - 2 do (
            m = M_mi;
            mShape = {numgens(target(m)), numgens(source(m))};
            mValences = sort(sumlist(entries(m), Axis => "row"));

            for ni from mi + 1 to #M - 1 do(
                if #positions(toSave, jj -> jj == ni) < 1 then (
                    continue;
                )
                else (
                    n = M_ni;
                    if n == m or n == -1*m then (
                        toSave = delete(ni, toSave);
                    )
                    else (
                        nShape = {numgens(target(m)), numgens(source(m))};
                        nValences = sort(sumlist(entries(n), Axis => "row"));
                        if nShape == mShape and nValences == mValences then (
                            possiblePermutations = combinations((0..(nShape#1-1)), nShape#1, R=>false);
                            for p in possiblePermutations do (
                                nPoss = n_p;
                                if isPerm(nPoss, m) or isPerm(nPoss, -1*m) then(
                                    toSave = delete(ni, toSave);
                                );
                            );
                        );
                    );
                );
            );
        );
        M_toSave
    )
    else M
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
