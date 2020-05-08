newPackage(
    "ToricQuiver",
    Headline => "Construction of flow polytopes and their associated quivers in arbitrary (up to computing power) dimension",
    Version => "0.0",
    Date => "November 20, 2019",
    Authors => {
        {Name => "Mary Barker",
         Email => "marybarker@wustl.edu",
         HomePage => "https://github.com/marybarker"}, 
        {Name => "Patricio Gallardo"}
    },
    PackageImports => {"Graphs"}
)
export {
-- Methods/Functions
    "sampleQuiver",
    "toricQuivers",
    "isTight",
    "subquivers",
    "subsetsClosedUnderArrows",
    "isStable",
    "isMaximal",
    "isAcyclic",
    "maximalUnstableSubquivers",
    "theta",
    "neighborliness",
    "flowPolytope",
    "graphEdges",
    "graphFromEdges",
    "wallType",
    "walls",
    "mergeOnVertex",
    "mergeOnArrow",
-- Options
    "Flow",
    "MatrixType",
    "AddOrientation",
    "Axis",
    "SavePath",
    "MaxSum",
    "MinSum",
    "Oriented",
    "RavelLoops",
    "Replacement",
    "EdgesAdded",
---- Quiver objects 
    "ToricQuiver",
    "toricQuiver"
}
protect Q0
protect Q1
protect flow
protect weights
protect connectivityMatrix

ToricQuiver = new Type of HashTable
toricQuiver = method(Options => {Flow=>"Canonical",MatrixType=>"Connectivity",AddOrientation=>false})

FlowCeil := 100;

-- construct ToricQuiver from connectivity matrix
toricQuiver(Matrix) := opts -> Q -> (
    F := 0.5*sumList(for x in entries(Q) list(for y in x list(abs(y))), Axis=>"Col");
    if opts.MatrixType == "Adjacency" then (
        Q = adjacencyToConnectivity(Q);
    );
    if opts.Flow == "Canonical" then (
        F = numColumns(Q):1;
    ) else if opts.Flow == "Random" then (
        F = for i in (0..#F - 1) list(random(FlowCeil));
        Q = Q*diagonalMatrix F;
    );
    new ToricQuiver from hashTable{
        connectivityMatrix=>Q,
        Q0=>toList(0..numRows(Q) - 1),
        Q1=>graphEdges(Q),
        flow=>F,
        weights=>sumList(entries(Q), Axis=>"Row")
    }
)

-- construct ToricQuiver from connectivity matrix and a flow
toricQuiver(Matrix, List) := opts -> (Q, F) -> (
    if opts.MatrixType == "Adjacency" then (
        Q = adjacencyToConnectivity(Q);
    );
    new ToricQuiver from hashTable{
        connectivityMatrix=>Q*diagonalMatrix F,
        Q0=>toList(0..numRows(Q) - 1),
        Q1=>graphEdges(Q),
        flow=>asList(F),
        weights=>sumList(entries(Q), Axis=>"Row")
    }
)

-- construct ToricQuiver from list of edges
toricQuiver(List) := E -> (
    Q := graphFromEdges(E, Oriented=>true);
    new ToricQuiver from hashTable{
        connectivityMatrix=>Q,
        Q0=>asList(0..numRows(Q) - 1),
        Q1=>E,
        flow=>asList(#E:1),
        weights=>sumList(entries(Q), Axis=>"Row")
    }
)

-- construct ToricQuiver from list of edges and a flow
toricQuiver(List, List) := (E, F) -> (
    Q := graphFromEdges(E, Oriented=>true)*diagonalMatrix F;
    new ToricQuiver from hashTable{
        connectivityMatrix=>Q,
        Q0=>toList(0..numRows(Q) - 1),
        Q1=>E,
        flow=>F,
        weights=>sumList(entries(Q), Axis=>"Row")
    }
)
-- subquiver of a ToricQuiver by taking a subset of the arrows, together with the vertices they adjoin
ToricQuiver _ List := (TQ, L) -> (
    M := matrix(for x in entries(TQ.connectivityMatrix_L) list(if sumList(x) != 0 then (x) else (continue;)));
    toricQuiver(M)
)
------------------------------------------------------------


------------------------------------------------------------
adjacencyToConnectivity = (A) -> (
    E := for i in (0..numRows(A) - 1) list(for j in (0..numColumns(A) - 1) list(if A_{j}^{i} != 0 then (i, j)));
    matrix(graphFromEdges(E), Oriented => true)

)
------------------------------------------------------------


------------------------------------------------------------
asList = x -> (
    if instance(x, List) then(
        return x
    )
    else if instance(x, Sequence) then(
        return toList(x)
    )
    else if instance(x, Set) then(
        return toList(x)
    )
    else
        return {x}
)
------------------------------------------------------------


------------------------------------------------------------
-- add all elements of a list x together, and specify Axis (row/col) if x is actually a matrix or list of lists -- 
sumList = {Axis => "None"} >> opts -> x -> (
    s := 0;
    if opts.Axis == "Row" then (
        s = flatten(for i in x list(sumList(i)));
    )
    else if opts.Axis == "Col" then (
       pivoted := entries(transpose(matrix(x)));
       s = flatten(for i in pivoted list(sumList(i)));
    )
    else (
        s = sum(asList(x));
    );
    return s
)
------------------------------------------------------------


------------------------------------------------------------
-- take all possible combinations of length k from list l -- 
-- optional arguments: 
-- -- Replacement(true/false) = with replacement
-- -- MinSum(numeric value) = exclude all combinations with sum below MinSum
-- -- MaxSum(numeric value) = exclude all combinations with sum above MaxSum
-- -- Order(true/false) = whether or not the ordering of combination values matters
combinations = {Replacement => true, MinSum => -1000, MaxSum => -1000, Order => true} >> opts -> (k, l) -> (
    combs := {};
    combs1 := {};
    combs2 := {};
    if k > 1 then (
        -- if we are using combinations with replacement -- 
        if opts.Replacement then (
           combs = flatten(join(for i in l list(for j from 0 to k - 1 list(i))));
           combs1 = unique(subsets(combs, k));
           combs2 = unique(subsets(combs, k));
           for i in combs2 do (combs1 = append(combs1, reverse(i)));
           combs = unique(combs1);
        )
        else (
           combs = flatten(for i in l list(i));
           combs1 = unique(subsets(combs, k));
           combs2 = unique(subsets(combs, k));
           for i in combs2 do (combs1 = append(combs1, reverse(i)));
           combs = unique(combs1);
        );
    )
    else
        combs = for i in l list(asList(i));

    -- if we are using restricting either by a minimum or maximum sum -- 
    if opts.MinSum != -1000 then (
       combs = for i in combs list(if sumList(i) < opts.MinSum then (continue;) else (i));
    );
    if opts.MaxSum != -1000 then (
       combs = for i in combs list(if sumList(i) > opts.MaxSum then (continue;) else (i));
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
    sortedVals := unique(sort(l));
    flatten(for i in sortedVals list(positions(l, x -> x == i)))
)
------------------------------------------------------------


------------------------------------------------------------
replaceInList = (i, v, l) -> (
    insert(i, v, drop(l, {i,i}))
)
------------------------------------------------------------


------------------------------------------------------------
isPermutation = (x, y) -> (
    toRet := false;

    xrows := entries(x);
    yrows := sort(entries(y));
    xcols := entries(transpose(x));
    ycols := entries(transpose(y));

    if #xrows == #yrows and #xcols == #ycols then (
        rs := toList(0..#xrows - 1);
        -- rowPermutations = permutations(rs);
        for rPerm in permutations(rs) do (
            if xrows_rPerm == yrows then (
                toRet = true;
                break;
            )
            else (
                xrowsP := matrix(xrows_rPerm);
                cs := toList(0..#xcols - 1);
                -- colPermutations := permutations(cs);
                for cPerm in permutations(cs) do (
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
-- get unique entries from list of undirected graphs -- 
unorientedUniqueUpToPermutation = x -> (
    if #x > 1 then (
        toSave := (0..(#x - 1));
        for i from 0 to #x - 2 do (
            for j from i + 1 to #x - 1 do (
                if #positions(toSave, a -> a == j) > 0 then (
                    if isPermutation(x#i, x#j) == false then (
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
allPossibleBaseGraphsForPair = x -> (
   g0 := x#0;
   g1 := x#1;

   -- get all possible columns for connectivity matrix (entries must be 0,1, or 2)
   possibleCols := combinations(g0, {0,1,2}, Replacement => true, MinSum=>2, MaxSum=>2);

   -- all combinations of columns to create rows
   rowCombs := combinations(g1, (0..(#possibleCols - 1)), Replacement=>true);
   candidateMats := for i in rowCombs list(for j in i list(asList(possibleCols#j)));
   candidateMats = for i in candidateMats list(transpose(matrix(i)));
   candidateMats = for i in candidateMats list(if min(sumList(entries(i), Axis=>"Row")) >= 3 then (i) else continue);
   candidateMats = unorientedUniqueUpToPermutation(candidateMats);
   return for m in candidateMats list (if isGraphConnected(m) then (m) else (continue;));
)
------------------------------------------------------------


------------------------------------------------------------
-- create all possible undirected base graphs for quiver generation in dimension d-- 
undirectedGraphs = (d) -> (
   graphPairs := apply((1..2*(d - 1)), i -> {i, i+d-1});
   connectivityMatrices := allPossibleBaseGraphsForPair \ graphPairs;
   flatten(connectivityMatrices)
)
------------------------------------------------------------


------------------------------------------------------------
-- yield the edges of a graph in the form of a list of pairs 
-- (v1, v2), where edge E is from v1 to v2
graphEdges = {Oriented => false, RavelLoops => false} >> opts -> (G) -> (
    E := {};
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
    lens := sortedIndices(for e in E list(-#e));
    return E_lens
)
------------------------------------------------------------


------------------------------------------------------------
-- yield the matrix rep of graph, given a list of edges as ordered 
-- pairs (this is the opposite of graphEdges() function. 
graphFromEdges = {Oriented => false} >> opts -> E -> (
    -- first, if oriented graph, then make sure this is reflected. 
    tailVal := 1;
    if opts.Oriented == true then (
        tailVal = -1;
    );

    nVerts := max(flatten(E))+1;
    cols := for i in E list(
        row := (nVerts:0);
        asList(replaceInList(i#0, tailVal, replaceInList(i#1, 1, row)))
    );
    transpose(matrix(cols))
)
------------------------------------------------------------


------------------------------------------------------------
edgesOutOfPoint = {Oriented => false} >> opts -> (p, E) -> (
    if opts.Oriented then (
        for i from 0 to #E - 1 list(e := E#i; if p != e#0 then (continue;) else (i, e))
    )
    else (
        for i from 0 to #E - 1 list(e := E#i; if (#positions(e, j -> j == p) < 1) then (continue;) else (i, e))
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- check if graph is connected
isGraphConnected = G -> (
    gEdges := graphEdges(G);
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
-- DFS search to find cycle in directed graph:
findCycleDFS = (startV, visited, E) -> (
    retVal := false;
    edgesOut := edgesOutOfPoint(startV, E, Oriented => true);
    for edge in edgesOut do (
        currentVisited := asList(visited);
        edgeVerts := edge#1;
        endV := edgeVerts#1;
        if visited#endV == 1 then (
            retVal = true;
            break;
        );
        currentVisited = replaceInList(endV, 1, visited);
        retVal = findCycleDFS(endV, currentVisited, E);
    );
    retVal
)
------------------------------------------------------------


------------------------------------------------------------
-- check if there exists a cycle in a (possibly unconnected)
-- oriented graph, passed in matrix form. 
existsOrientedCycle = (G) -> (
    retVal := false;
    E := graphEdges(G, Oriented => true);
    V := asList(0..numRows(G)-1);
    for firstV in V do (
        visited := insert(firstV, 1, asList((#V - 1):0));
        result := findCycleDFS(firstV, visited, E);
        if result then (
            retVal = true;
            break;
        )
    );
    retVal
)
------------------------------------------------------------


------------------------------------------------------------
isAcyclic = method()
isAcyclic(Matrix) := Q -> (
    not existsOrientedCycle(Q)
)
isAcyclic(ToricQuiver) := Q -> (
    not existsOrientedCycle(Q.connectivityMatrix)
)
------------------------------------------------------------


------------------------------------------------------------
-- check if there exists a path between p and q by appending 
-- edges in E(which is a list of pairs (v1, v2). 
-- optional arguments: 
-- -- Oriented(true/false) = whether or not the graph should be oriented
-- -- SavePath(true/false) = whether or not to return the edges involved in the path
-- -- EdgesAdded(list) = internal mechanism for computing for SavePath

------------------------------------------------------------
isPathBetween = {Oriented => false, SavePath => false, EdgesAdded => {}} >> opts -> (p, q, E) -> (
    ifPath := false;
    existsPath := false;
    currentEdges := {};
    pathsToSee := edgesOutOfPoint(p, E, Oriented => opts.Oriented);

    for edge in pathsToSee do (
        --- get the edge index and enpoints
        i := edge#0;
        e := edge#1;
        v := e#1;
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
            thisPath := {};
            remainingEdges := for j from 0 to #E - 1 list(if j == i then (continue;) else E#j);

            if opts.SavePath then (
                (ifPath, thisPath) = isPathBetween(v, q, remainingEdges, Oriented => opts.Oriented, SavePath => true, EdgesAdded => currentEdges);
            )
            else (
                ifPath = isPathBetween(v, q, remainingEdges, Oriented => opts.Oriented, EdgesAdded => currentEdges);
            );
            if ifPath then (
                existsPath = true;
                currentEdges = currentEdges | thisPath;
                break;
            );
        );
    );
    if opts.SavePath then (
        return (existsPath, currentEdges);
    )
    else (
        return existsPath;
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- checks if there is a cycle containing given edge. 
isEdgeInCycle = (i, E) -> (
    if #E > 1 then (
        e := E#i;
        if #e > 1 then (
            p := e#0;
            q := e#1;
        )
        else (
            p = e#0;
            q = e#0;
        );
        indicesToSave := drop(toList(0..(#E-1)), {i,i});
        isPathBetween(p, q, E_indicesToSave)
    )
    else (
        false
    )
)
------------------------------------------------------------


------------------------------------------------------------
splitLoops = m -> (
    Es := graphEdges(m);
    nVerts := numRows(m);
    loopsBroken := for i from 0 to #Es - 1 list(
        e := Es#i;
        if (#e < 2) or #(delete(e#0, e)) < 1 then (
            i
        )
        else (continue;)
    );
    altLB := flatten(for i in loopsBroken list(i));
    newEdges := flatten(for i from 0 to #Es - 1 list (
        e := Es#i;
        if #loopsBroken != #delete(i, loopsBroken) then (
            p := position(loopsBroken, x -> x == i);
            altLB = append(replaceInList(p, loopsBroken_p + p, altLB), loopsBroken_p + p + 1);
            {{e#0, nVerts+p}, {nVerts+p, e#0}}
        )
        else (
            {e}
        )
    ));
    (graphFromEdges(newEdges), altLB)
)
------------------------------------------------------------


------------------------------------------------------------
splitEdges = (m, E) -> (
    Es := graphEdges(m);
    nVerts := numRows(m);

    for i from 0 to #E - 1 do (
        ei := E#i;
        e := Es#ei;
        Es = append(replace(i, {e#0, nVerts+i}, Es), {nVerts+i, e#1});
    );
    graphFromEdges(Es)
)
------------------------------------------------------------


------------------------------------------------------------
-- generate the undirected graphs with admissable #vertices and #arrows to generate quivers
undirectedBaseGraphs = (n) -> (
    for m in undirectedGraphs(n) list(
        es := graphEdges(m, RavelLoops => true);
        if #es > 0 then (
            inCycle := for i from 0 to #es - 1 list( isEdgeInCycle(i, es) );
            if all(inCycle, i -> i == true) then (m) else (continue;)
        )
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- remove all loops from a given graph given in matrix form
-- and return the set of all loop-free copies with every 
-- combination of the remaining edges split
splitLoopsAndEdges = m -> (
    origNumEdges := #graphEdges(m);
    A := splitLoops(m);
    M := A#0;
    loops := A#1;
    originalLooplessEdges := set(0..origNumEdges - 1) - set(loops);


    if #originalLooplessEdges > 0 then (
        otherEdgeCombs := for i from 1 to #originalLooplessEdges list(
            asList(combinations(i, asList(originalLooplessEdges), Replacement => false, Order => false))
        );
        outputs := flatten(for i from 0 to #originalLooplessEdges - 1 list(
            for comb in otherEdgeCombs#i list(
                splitEdges(M, asList(comb))
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
allPossibleOrientations = m -> (
    rows := entries(m);

    -- first make all vertices of valence 2 into sinks. 
    mat := matrix(for i in rows list(
        if sumList(i) == 2 then (-1*i) else (i)
    ));
    cols := entries(transpose(mat));
    es := graphEdges(mat);

    -- now take every combination of +-1 for the columns that do 
    -- not yet add up to 0
    columnsToChange := for i from 0 to #cols - 1 list(
        if sumList(cols#i) == 0 then (continue;) else (i)
    );

    if #columnsToChange > 0 then (
        columnChoices := for i in columnsToChange list(
            edge := es#i;
            tail := edge#0;
            head := edge#1;
            {replaceInList(tail, -1, cols#i), replaceInList(head, -1, cols#i)}
        );

        combinationsOfChoices := combinations(#columnsToChange, {0, 1}, Replacement=>true, Order=>false);
        for choice in combinationsOfChoices list(-1*transpose(matrix(
            for c from 0 to #cols - 1 list(
                if #columnsToChange == #delete(c, columnsToChange) then (
                    cols#c
                )
                else (
                    p := position(columnsToChange, i -> i == c);
                    colOpts := columnChoices#p;
                    choiceForCol := choice#p;
                    colOpts#choiceForCol
                )
            )
        )))
    )
    else (
        {-1*mat}
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- take a list of quivers(directed graphs) and return the ones
-- that do not contain any oriented cycles. 
removeOrientedCycles = (l) -> (
    for m in l list(
        if existsOrientedCycle(m) then (
            continue;
        )
        else (
            m
        )
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- take a list of quivers(directed graphs) and return the 
-- unique up to isomorphism list 
uniqueUpToQuiverIsomorphism = M -> (
    if #M > 1 then (
        toSave := asList(0..(#M - 1));
        for mi from 0 to #M - 2 do (
            m := M_mi;
            mShape := {numgens(target(m)), numgens(source(m))};
            mValences := sort(sumList(entries(m), Axis => "Row"));

            for ni from mi + 1 to #M - 1 do(
                if #positions(toSave, jj -> jj == ni) < 1 then (
                    continue;
                )
                else (
                    n := M_ni;
                    if n == m or n == -1*m then (
                        toSave = delete(ni, toSave);
                    )
                    else (
                        nShape := {numgens(target(m)), numgens(source(m))};
                        nValences := sort(sumList(entries(n), Axis => "Row"));
                        if nShape == mShape and nValences == mValences then (
                            possiblePermutations := combinations(nShape#1, (0..(nShape#1-1)), Replacement=>false);
                            for p in possiblePermutations do (
                                nPoss := n_p;
                                if isPermutation(nPoss, m) or isPermutation(nPoss, -1*m) then(
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
-- generates the Toric Quivers in dimension n
toricQuivers = n -> (
    gs := undirectedBaseGraphs(n);
    gs = flatten(for i in gs list(splitLoopsAndEdges(i)));
    gs = flatten(for i in gs list(allPossibleOrientations(i)));
    qs := removeOrientedCycles(gs);
    qs = uniqueUpToQuiverIsomorphism(qs);
    for q in qs list(
        toricQuiver(q)
    )
)
------------------------------------------------------------


------------------------------------------------------------
getFirstGraphFromPair = (g0, g1) -> (
    -- get all possible columns for connectivity matrix (entries must be 0,1, or 2)
    possibleCols := combinations(g0, {0,1,2}, Replacement => true, MinSum=>2, MaxSum=>2);

    -- all combinations of columns to create rows
    rowCombs := combinations(g1, (0..(#possibleCols - 1)), Replacement=>true);
    M := 0;
    firstG := 0;
    for i in rowCombs do(
        M = transpose(matrix(for j in i list(asList(possibleCols#j))));
        if min(sumList(entries(M), Axis=>"Row")) >= 3 then (
            firstG = M; 
            break;
        ) else (continue;)
    );
    firstG
)
------------------------------------------------------------


------------------------------------------------------------
-- extract a sample quiver in dimension n 
sampleQuiver = {Flow => "Canonical"} >> opts -> (n) -> (
    metSq := false;
    sq := {};
    tries := 1;
    g0 := 0;
    g1 := 0;
    possibleCols := 0;
    rowCombs := 0;
    M := 0;

    while (not metSq) and (tries < 2*n - 1) do (
        g0 = tries;
        g1 = g0 + n - 1;
        tries = tries + 1;
        possibleCols = combinations(g0, {0,1,2}, Replacement=>true, MinSum=>2, MaxSum=>2);
        rowCombs = combinations(g1, (0..(#possibleCols - 1)), Replacement=>true);

        for rc in rowCombs do(
            if metSq then break;
            M = transpose(matrix(for j in rc list(asList(possibleCols#j))));
            if min(sumList(entries(M), Axis=>"Row")) >= 3 then (
                for i in splitLoopsAndEdges(M) do (
                    if metSq then break;
                    for j in allPossibleOrientations(i) do (
                        if (not existsOrientedCycle(j)) then (
                            sq = toricQuiver(matrix(j));
                            metSq = true;
                            break;
                        )
                    )
                )
            )
        )
    );
    if opts.Flow == "Canonical" then (
        sq
    ) else (
        toricQuiver(sq.connectivityMatrix, Flow=>"Random")
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- yield the subquivers of a given quiver Q
subquivers = method(Options => {Format => "quiver"})
subquivers Matrix := opts -> Q -> (
    numArrows := numColumns(Q);
    arrows := 0..(numArrows - 1);

    flatten(
        for i from 1 to numArrows - 1 list (
            for c in combinations(i, arrows, Order => false, Replacement => false) list (
                if opts.Format == "indices" then (
                    c
                ) else (
                    toricQuiver(Q)_c
                )
            )
        )
    )
)
subquivers ToricQuiver := opts -> Q -> (
    numArrows := #Q.Q1;
    arrows := Q.Q1;

    flatten(
        for i from 1 to numArrows - 1 list (
            for c in combinations(i, arrows, Order => false, Replacement => false) list (
                if opts.Format == "indices" then (
                    c
                ) else (
                    Q_c
                )
            )
        )
    )
)
------------------------------------------------------------


------------------------------------------------------------
-- list the subsets of a quiver Q that are closed under arrows
subsetsClosedUnderArrows = method()
subsetsClosedUnderArrows Matrix := (Q) -> (
    numVertices := numRows(Q);
    currentVertices := 0..(numVertices - 1);
    Qt := transpose(Q);

    flatten(for i from 1 to numVertices - 1 list(
        for c in combinations(i, currentVertices, Order => false, Replacement => false) list(
            sQ := entries(Qt_c);
            if all(sumList(sQ, Axis => "Row"), x -> x >= 0) then (
                c
            )
            else(
                continue;
            )
        )
    ))
)
subsetsClosedUnderArrows ToricQuiver := (Q) -> (
    subsetsClosedUnderArrows(Q.connectivityMatrix)
)
------------------------------------------------------------


------------------------------------------------------------
-- return ordered list of the weights for the vertices of quiver Q
theta = method()
theta(ToricQuiver) := Q -> (
    Q.weights
)
theta(Matrix) := Q -> (
    sumList(entries(Q), Axis => "Row")
)
------------------------------------------------------------


------------------------------------------------------------
isStable = method()
isStable(Matrix, List) := (Q, subQ) -> (
    -- get the vertices in the subquiver
    subQVertices := positions(entries(Q_subQ), x -> any(x, y -> y != 0));
    -- weights of the original quiver
    Qtheta := theta(Q);
    -- inherited weights on the subquiver
    weights := Qtheta_subQVertices;
    subMat := Q_subQ;
    tSubMat := transpose(subMat);
    subMat = transpose(tSubMat_subQVertices);

    sums := asList(
        for subset in subsetsClosedUnderArrows(subMat) list(
            sumList(weights_subset)
        )
    );
    all(sums, x -> x > 0)
)
isStable(ToricQuiver, List) := (Q, subQ) -> (
    isStable(Q.connectivityMatrix, subQ)
)
------------------------------------------------------------


------------------------------------------------------------
unstableSubquivers = method()
unstableSubquivers Matrix := Q -> (
    numArrows := numColumns(Q);
    arrows := asList(0..numArrows - 1);

    L := flatten(for i from 1 to numArrows - 1 list (
        combinations(i, arrows, Replacement => false) 
    ));

    for sQ in L list(
        if not isStable(Q, asList(sQ)) then (
            sQ
        ) else (
            continue;
        )
    )
)
unstableSubquivers ToricQuiver := Q -> (
    unstableSubquivers(Q.connectivityMatrix)
)
------------------------------------------------------------


------------------------------------------------------------
isProperSubset = (Q1, Q2) -> (
    if set(Q1) === set(Q2) then (
        false
    ) else (
        isSubset(set(Q1), set(Q2))
    )
)
------------------------------------------------------------


------------------------------------------------------------
isMaximal = method()
isMaximal(Matrix, List) := (Q, Qlist) -> (
    returnVal := true;
    for Q2 in Qlist do (
        if isProperSubset(Q, Q2) then (
            returnVal = false;
        );
    );
    returnVal
)
isMaximal(ToricQuiver, List) := (Q, Qlist) -> (
    isMaximal(Q.connectivityMatrix, Qlist)
)
------------------------------------------------------------


------------------------------------------------------------
maximalUnstableSubquivers = (Q) -> (
    unstableList := unstableSubquivers(Q);
    for subQ1 in unstableList list (
        isMaximal := true;
        for subQ2 in unstableList do (
            if isProperSubset(subQ1, subQ2) then (
                isMaximal = false;
            );
        );
        if isMaximal then (
            subQ1
        ) else (
            continue;
        )
    )
)
------------------------------------------------------------


------------------------------------------------------------
isTight = method()
isTight(Matrix) := Q -> (
    numArrows := numColumns(Q);
    all(maximalUnstableSubquivers(Q), x -> #x != (numArrows - 1))
)

isTight(ToricQuiver) := Q -> (
    numArrows := #Q#Q1;
    all(maximalUnstableSubquivers(Q.connectivityMatrix), x -> #x != (numArrows - 1))
)
------------------------------------------------------------


------------------------------------------------------------
neighborliness = method()
neighborliness Matrix := (Q) -> (
    numArrows := numColumns(Q);
    maxUnstables := maximalUnstableSubquivers(Q);

    k := max(
        for sQ in maxUnstables list(
            numArrows - #sQ
        )
    );
    k
)
neighborliness ToricQuiver := (Q) -> (
    neighborliness(Q.connectivityMatrix)
)
------------------------------------------------------------


------------------------------------------------------------
wallType = method()
wallType(Matrix, List) := (Q, Qp) -> (
    tp := sum(for x in sumList(Q^Qp, Axis=>"Col") list(if x < 0 then (1) else (continue;)));
    tm := sum(for x in sumList(Q^Qp, Axis=>"Col") list(if x > 0 then (1) else (continue;)));
    (tp, tm)
)
wallType(ToricQuiver, List) := (Q, Qp) -> (
    wallType(Q.connectivityMatrix, Qp)
)
------------------------------------------------------------


------------------------------------------------------------
walls = method()
walls(Matrix) := (Q) -> (
    nv := numRows(Q);
    nvSet := set(0..nv - 1);
    subs := (1..ceiling(nv/2));

    Qms := flatten(for i from 1 to ceiling(nv/2) list (
        combinations(i, asList(nvSet), Replacement => false)
    ));

    Qedges := graphEdges(Q, Oriented=>true);
    for Qm in Qms list(
        mSums := sumList(Q^Qm, Axis=>"Col");
        QmEdgeIndices := for s in (0..#mSums - 1) list(if (mSums_s == 0) then (s) else (continue;));

        if (#Qm < 2) or (isGraphConnected(Q^Qm_QmEdgeIndices)) then (
            Qp := asList(nvSet - set(Qm));
            pSums := sumList(Q^Qp, Axis=>"Col");
            QpEdgeIndices := for s in (0..#pSums - 1) list(if (pSums_s == 0) then (s) else (continue;));
            if (#Qp < 2) or (isGraphConnected(Q^Qp_QpEdgeIndices)) then (
               (Qp, wallType(Q, Qp))
            )
        )
    )
)
walls(ToricQuiver) := (Q) -> (
    walls(Q.connectivityMatrix)
)
------------------------------------------------------------


------------------------------------------------------------
-- Returns a spanning tree(the first one that is encountered) of 
-- the quiver Q with |Q_1| - |Q_0| + 1 edges removed. 
-- NOTE: if such a spanning tree is not possible, then it returns empty lists
--
-- input: 
--     - Q: Matrix representation of quiver
-- outputs:
--     - edges_kept(list of tuples): list of the edges in the spanning tree
--     - edges_removed(list of tuples)): list of the edges in the complement of the spanning tree
--
spanningTree = (Q) -> (
    Q0 := numRows(Q);
    Q1 := numColumns(Q);

    --  edges of quiver Q represented as a list of tuples
    allEdges := graphEdges(Q, Oriented => true);
    allNodes := asList(0..Q0-1);

    -- number of edges to remove from spanning tree
    d := Q1 - Q0 + 1;

    dTuplesToRemove := combinations(d, asList(0..#allEdges-1), Replacement => false);
    for dTuple in dTuplesToRemove do (
        edgesKept := drop(allEdges, dTuple);
        edgesRemoved := allEdges_dTuple;

        reducedG := transpose(matrix(for e in edgesKept list(
            t := e#0;
            h := e#1;
            localE := asList(Q0:0);
            localE = replaceInList(h,  1, localE);
            localE = replaceInList(t, -1, localE);
            localE
        )));
        if isGraphConnected(reducedG) then (
            return (edgesKept, edgesRemoved);
        );
    );
    return ({}, {});
)
------------------------------------------------------------


------------------------------------------------------------
isIn = (v, l) -> (
    p := positions(l, x -> x == v);
    #p > 0
)
------------------------------------------------------------


------------------------------------------------------------
-- gives the edges that comprise an undirected cycle in the graph G, 
-- (which is assumed to contain a single cycle) and returns the ordered cycle

--  input: G(list of tuples): edges of graph G
--  output: cycle(list of tuples): tuple representation of the edges contained in the cycle
primalUndirectedCycle = (G) -> (
    for i from 0 to #G - 1 do (
        edge := G#i;
        (isCycle, cycle) := isPathBetween(edge#1, edge#0, drop(G, {i, i}), 
                                          Oriented => false, SavePath => true, EdgesAdded => {edge});
        if isCycle then (
            edgeIndices := {};
            metEdges := {};

            for cE in cycle do (
                for gI in asList(0..#G - 1) do (
                    if isIn(gI, metEdges) then (
                        continue;
                    ) else (
                        gE := G#gI;
                        if (gE#0 == cE#0) and (gE#1 == cE#1) then (
                            metEdges = metEdges | {gI};
                            edgeIndices = edgeIndices | {gI};
                        )
                        else if (gE#1 == cE#0) and (gE#0 == cE#1) then (
                            metEdges = metEdges | {gI};
                            edgeIndices = edgeIndices | {-(gI+1)};
                        );
                    );
                );
            );
        );
        return edgeIndices;
    );
    return {};
)
------------------------------------------------------------


------------------------------------------------------------
flowPolytope = method()
flowPolytope(Matrix) := (Q) -> (

    (sT, removedEdges) := spanningTree(Q);
    es := sT | removedEdges;
    Ws := #es:1;

    f := for i from 0 to #removedEdges - 1 list(
        edge := removedEdges#i;
        cycle := primalUndirectedCycle(sT | {edge});

        cycle = for x in cycle list(
            if x == #sT then (x+i) else if x == -(#sT + 1) then (-#sT - i - 1) else (x)
        );

        fi := #es:0;
        for j in cycle do (
            if j >= 0 then (
                fi = replaceInList(j, Ws#j, fi)
            ) else (
                k := -(1 + j);
                fi = replaceInList(k, -Ws#k, fi)
            );
        );
        fi
    );
    for j from 0 to #es - 1 list(
        for i from 0 to #removedEdges - 1 list(
            ff := f#i;
            ff#j
        )
    )
)
flowPolytope(ToricQuiver) := (Q) -> (
    flowPolytope(Q.connectivityMatrix)
)
------------------------------------------------------------


------------------------------------------------------------
mergeOnVertex = method()
mergeOnVertex(Matrix, ZZ, Matrix, ZZ) := (Q1, v1, Q2, v2) -> (
    nrow := numRows(Q1) + numRows(Q2) - 1;
    ncol := numColumns(Q1) + numColumns(Q2);

    i1 := asList(join(drop(0..numRows(Q1) - 1, {v1, v1}), {v1}));
    i2 := asList(join({v2}, drop(0..numRows(Q2) - 1, {v2, v2})));
    Q1 = transpose(Q1^i1);
    Q2 = transpose(Q2^i2);

    paddingSize := 0;
    matrix(
        for row from 0 to nrow - 1 list(
            if row < (numColumns(Q1) - 1) then (
                paddingSize = ncol - numColumns(Q1) - 1;
                join(entries(Q1)_row, paddingSize:0)
            ) else if row < numColumns(Q1) then (
                join(entries(Q1)_row, entries(Q2)_0)
            ) else (
                j := row - numColumns(Q1) + 1;
                paddingSize = ncol - numRows(Q2);
                asList(join(paddingSize:0, entries(Q2)_j))
            )
        )
    )
)
mergeOnVertex(ToricQuiver, ZZ, Matrix, ZZ) := (Q1, v1, Q2, v2) -> (
    mergeOnVertex(Q1.connectivityMatrix, v1, Q2, v2)
)
mergeOnVertex(Matrix, ZZ, ToricQuiver, ZZ) := (Q1, v1, Q2, v2) -> (
    mergeOnVertex(Q1, v1, Q2.connectivityMatrix, v2)
)
mergeOnVertex(ToricQuiver, ZZ, ToricQuiver, ZZ) := (Q1, v1, Q2, v2) -> (
    mergeOnVertex(Q1.connectivityMatrix, v1, Q2.connectivityMatrix, v2)
)
------------------------------------------------------------


------------------------------------------------------------
mergeOnArrow = method()
mergeOnArrow(Matrix, ZZ, Matrix, ZZ) := (Q1, a1, Q2, a2) -> (
    nrow := numRows(Q1) + numRows(Q2) - 2;
    ncol := numColumns(Q1) + numColumns(Q2) - 1;

    q1E := asList(graphEdges(Q1, Oriented=>true))_a1;
    q2E := asList(graphEdges(Q2, Oriented=>true))_a2;

    c1 := asList(join(drop(0..numColumns(Q1) - 1, {a1, a1}), {a1}));
    c2 := asList(drop(0..numColumns(Q2) - 1, {a1,a1}));

    r1 := asList(join(asList(set(0..numRows(Q1) - 1) - set(q1E)), q1E));
    r2 := asList(join(q2E, asList(set(0..numRows(Q2) - 1) - set(q2E))));

    Q1 = transpose((Q1^r1)_c1);
    Q2 = transpose((Q2^r2)_c2);

    paddingSize := 0;
    matrix(
        for row from 0 to nrow - 1 list(
            if row < (numColumns(Q1) - 2) then (
                paddingSize = ncol - numColumns(Q1) - 1;
                join(entries(Q1)_row, asList(paddingSize:0))

            ) else if row < numColumns(Q1) then (
                join(entries(Q1)_row, entries(Q2)_(2 + row - numColumns(Q1)))
            ) else (
                j := row - numColumns(Q1) + 2;
                paddingSize = ncol - numRows(Q2);
                asList(join(paddingSize:0, entries(Q2)_j))
            )
        )
    )
)
mergeOnArrow(ToricQuiver, ZZ, Matrix, ZZ) := (Q1, a1, Q2, a2) -> (
    mergeOnArrow(Q1.connectivityMatrix, a1, Q2, a2)
)
mergeOnArrow(Matrix, ZZ, ToricQuiver, ZZ) := (Q1, a1, Q2, a2) -> (
    mergeOnArrow(Q1, a1, Q2.connectivityMatrix, a2)
)
mergeOnArrow(ToricQuiver, ZZ, ToricQuiver, ZZ) := (Q1, a1, Q2, a2) -> (
    mergeOnArrow(Q1.connectivityMatrix, a1, Q2.connectivityMatrix, a2)
)
------------------------------------------------------------
