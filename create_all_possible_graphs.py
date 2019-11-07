import graph_cal as gc
import numpy as np
import random
import os
from itertools import combinations
import networkx as nx
import matplotlib.pyplot as plt

n = int(input("Enter value for d -> "))

folder = "./outputs/d=%d/"%n
filename = folder+"step"

def read_step_file(filename):
    graphs = []
    with open(filename, "r") as f:
        A = []
        for line in f.readlines():
            line = line.rstrip('\n')
            if len(line) > 0:
                A.append([int(x) for x in line.split(' ')])
            else:
                if len(A) > 0:
                    graphs.append(np.matrix(A))
                    A = []
    return gc.unique_up_to_isomorphism(graphs)


if os.path.exists(filename+"1"):
    graphs = read_step_file(filename+"1")
else:
    graphs = gc.Step1(n)
print("step 1: ")
with open(filename+"1","w") as f:
    for g in graphs:
        for line in g:
            np.savetxt(f, line, fmt='%d')
        f.write("\n")




step2_graphs, loops_broken = zip(*[gc.Step2(A) for A in graphs])
print("step 2: ")
with open(filename+"2","w") as f:
    for g in step2_graphs:
        for line in g:
            np.savetxt(f, line, fmt='%d')
        f.write("\n")




#step3_graphs = [gc.Step3(A, edge_list) for edge_list in combinations(range(A.shape[1]), y) for y in range(A.shape[1]) for A in step2_graphs]
step3_graphs = []
for i, A in enumerate(step2_graphs):
    n_edges = A.shape[1]
    edges_to_break = [x for x in range(n_edges) if not (x in loops_broken[i])]
    if len(edges_to_break) > 0:
        for y in range(len(edges_to_break)+1):
            for edge_list in combinations(edges_to_break, y):
                step3_graphs.append(gc.Step3(A, list(edge_list)))
    else:
        step3_graphs.append(A)
step3_graphs = gc.unique_up_to_isomorphism(step3_graphs)

print("step 3: ")
with open(filename+"3","w") as f:
    for g in step3_graphs:
        for line in g:
            np.savetxt(f, line, fmt='%d')
        f.write("\n")




step4_graphs = []
for A in step3_graphs:
    step4_graphs.extend(gc.Step4(A))
step4_graphs = gc.unique_up_to_isomorphism(step4_graphs)
#step4_graphs = gc.unique_up_to_isomorphism([y for y in gc.Step4(A) for A in step3_graphs])
print("step 4: ")
with open(filename+"4","w") as f:
    for g in step4_graphs:
        for line in g:
            np.savetxt(f, line, fmt='%d')
        f.write("\n")




step5_graphs = gc.Step5(step4_graphs)
print("step 5: ")
with open(filename+"5","w") as f:
    for g in step5_graphs:
        for line in g:
            np.savetxt(f, line, fmt='%d')
        f.write("\n")


flow_polytopes = []
for i, g in enumerate(step5_graphs):
    flow_polytopes.append(gc.flow_polytope(g))
with open(filename+"_polytope","w") as f:
    for v in flow_polytopes:
        np.savetxt(f, v, fmt='%d')
        f.write('\n')

# now generate graphic that shows growth of 
G = nx.Graph()
G.add_nodes_from(["1_%d"%x for x in range(len(graphs))])

nums = np.zeros(len(graphs))
for n1, graph in enumerate(graphs):
    all_of_them = []
    g2s, loops_broken = gc.Step2(graph)

    for n2, g2 in enumerate([g2s]):
        G.add_edge("1_%d"%n1, "1_%d_2_%d"%(n1, n2), color="r")

        step3_graphs = []
        n_edges = g2.shape[1]
        edges_to_break = [x for x in range(n_edges) if not (x in loops_broken)]
        if len(edges_to_break) > 0:
                for y in range(len(edges_to_break)+1):
                    for edge_list in combinations(edges_to_break, y):
                        step3_graphs.append(gc.Step3(g2, list(edge_list)))
        else:
            step3_graphs.append(g2)
        g3s = gc.unique_up_to_isomorphism(step3_graphs)

        for n3, g3 in enumerate(g3s):
            G.add_edge("1_%d_2_%d"%(n1,n2), "1_%d_2_%d_3_%d"%(n1,n2,n3), color="g")
            g4s = gc.Step4(g3)

            for n4, g4 in enumerate(g4s):
                G.add_edge("1_%d_2_%d_3_%d"%(n1,n2,n3), "1_%d_2_%d_3_%d_4_%d"%(n1,n2,n3,n4), color="b")

                g5s = [x for x in g4s if not gc.exists_cycle(gc.edges_of_graph(x, True), range(x.shape[0]))]
                for n5, g5 in enumerate(g5s):
                    G.add_edge("1_%d_2_%d_3_%d_4_%d"%(n1,n2,n3,n4), "1_%d_2_%d_3_%d_4_%d_5_%d"%(n1,n2,n3,n4,n5), color="b")
                all_of_them.extend(g5s)
    nums[n1] = len(gc.unique_up_to_isomorphism(all_of_them))

plt.scatter(range(len(nums)), nums)
plt.savefig(folder+"graph_growth_rate.png")
plt.clf()

cs = [len(x.split('_')) for x in G.nodes]
nx.draw(G, node_size=20, node_color=cs)
plt.savefig(folder+"graph_of_connections.png")
