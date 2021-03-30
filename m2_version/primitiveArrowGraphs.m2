loadPackage "ThinSincereQuivers"
needsPackage "NormalToricVarieties";
needsPackage "Graphs";
needsPackage "Polyhedra"

shareAFacet = (C1, C2) -> (

    for f1 in facesAsCones(1, C1) do (
        for f2 in facesAsCones(1, C2) do (
            if (f1 == f2) then (
                return true;
            );
        );
    );
    return false;
)

chamberGraph = Q -> (
    print("toric quiver", Q);
    CQ := coneSystem Q;
    ths := referenceThetas CQ;
    fps := apply(ths, x-> transpose matrix flowPolytope(x, Q));
    nts := apply(fps, x-> normalToricVariety x);
    tif := apply(nts, x-> isFano x);
    nvs := apply(fps, x-> #latticePoints convexHull x);
    pgs := apply(nts, x-> picardGroup x);
    edges := {};

    for ic1 in (0..#CQ - 1) do(

        c1 := CQ#ic1;
        for ic2 in (ic1 + 1..#CQ - 1) do(

            c2 := CQ#ic2;
            if shareAFacet(c1, c2) then (
                edges = edges | {{ic1, ic2}};
            );
        );
    );
    verts := apply(0..#ths - 1, x-> {x, nvs#x, tif#x, fps#x, nvs#x, pgs#x});
    return {verts, edges, max(nvs) - min(nvs)};
)

ThreePrimitiveArrowQuivers := {
    --{{0,1},{0,2},{0,3}},
    --{{1,0},{2,0},{3,0}},
    --{{0,1},{1,2},{2,3}},
    --{{0,1},{1,2},{2,3},{0,2}},
    {{0,1},{1,2},{2,3},{0,3}},
    {{0,1},{1,2},{2,3},{1,3}},
    {{0,1},{1,2},{2,3},{0,2},{0,3}},
    {{0,1},{1,2},{2,3},{0,2},{1,3}},
    {{0,1},{1,2},{2,3},{0,3},{1,3}},
    {{0,1},{1,2},{2,3},{0,2},{0,3},{1,3}},
    --{{0,1},{1,2},{3,2}},
    --{{0,1},{1,2},{3,2},{0,2}},
    {{0,1},{1,2},{3,2},{0,3}},
    {{0,1},{1,2},{3,2},{0,2},{0,3}},
    --{{1,0},{1,2},{2,3}},
    --{{1,0},{1,2},{2,3},{1,3}},
    --{{1,0},{1,2},{2,3},{2,0}},
    {{1,0},{1,2},{2,3},{1,3},{2,0}},
    --{{0,1},{0,2},{3,0}},
    --{{0,1},{0,2},{3,0},{3,1}},
    --{{0,1},{0,2},{3,0},{3,2}},
    {{0,1},{0,2},{3,0},{3,1},{3,2}},
    --{{0,1},{2,0},{3,0},{3,1}},
    --{{0,1},{2,0},{3,0},{2,1}},
    {{0,1},{2,0},{3,0},{3,1},{2,1}}
};



fname = "threePrimitiveArrows";
for Q in ThreePrimitiveArrowQuivers do(
    TQ := toricQuiver(Q);
    CG := chamberGraph(TQ);
    fname << Q << endl;
    for v in CG#0 do(
        fname << v << endl;
        print(v);
    );
    for e in CG#1 do(
        fname << e << endl;
        print(e);
    );
    fname << CG#2 << endl;
    fname << endl;
);
fname << close;



fname = "lookAtAddingArrows";
Q = first ThreePrimitiveArrowQuivers;

for iq in (0..#Q - 1) do(
    i = Q#iq;

    TQ := toricQuiver(Q | {i});
    CG := chamberGraph(TQ);
    fname << quiverEdges TQ << endl;

    print(i, iq, TQ);

    for v in CG#0 do(
        fname << v << endl;
        print(v);
    );
    for e in CG#1 do(
        fname << e << endl;
        print(e);
    );
    fname << CG#2 << endl;
    fname << endl;

    for jq from iq to #Q - 1 do(
        j = Q#jq;

        TQ := toricQuiver(Q | {i} | {j});
        CG := chamberGraph(TQ);
        fname << quiverEdges TQ << endl;
        for v in CG#0 do(
            fname << v << endl;
            print(v);
        );
        for e in CG#1 do(
            fname << e << endl;
            print(e);
        );
        fname << CG#2 << endl;
	fname << endl;

        for kq from jq to #Q - 1 do(
            k = Q#kq;
            TQ := toricQuiver(Q | {i} | {j} | {k});
            CG := chamberGraph(TQ);
            fname << quiverEdges TQ << endl;
            for v in CG#0 do(
                fname << v << endl;
                print(v);
            );
            for e in CG#1 do(
                fname << e << endl;
                print(e);
            );
            fname << CG#2 << endl;
            fname << endl;
        );
    );
    fname << "next round" << endl << endl;
);
fname << close;


