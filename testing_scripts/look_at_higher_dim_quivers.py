import sys
sys.path.append("..")
import ToricQuiver as tq
import numpy as np

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
    return tq.unoriented_unique_up_to_isomorphism(graphs)

def write_graphs_to_file(gs, filename):
    with open(filename, "w") as f:
        for g in gs:
            f.write(g[1] + ' \n')
            for line in g[0]:
                np.savetxt(f, line, fmt='%d')
            f.write("\n")

def is_in_list(g, gs):
    rVal = False
    for g1 in gs:
        if np.all(g1 == g):
            rVal = True
            break
    return rVal


d1 = 1
d2 = 2

g1s = [np.matrix([[-1, -1], [1, 1]])]
#g1s = [tq.standard_form(x) for x in read_step_file("../outputs/d=%d/step5"%d1)]
g2s = [tq.standard_form(x) for x in read_step_file("../outputs/d=%d/step5"%d2)]
g3s = [tq.standard_form(x) for x in read_step_file("../outputs/d=%d/step5"%(d1+d2))]
unique_list = []

for i1, g1 in enumerate(g1s):
    g1_edges = tq.edges_of_graph(g1s[i1])

    for e1 in np.unique(g1_edges):
        for i2, g2 in enumerate(g2s):
            g2_edges = tq.edges_of_graph(g2)

            for e2 in np.unique(g2_edges):
                g3 = tq.standard_form(np.matrix(tq.merge_on_arrow(g1, e1, g2, e2)))

                if not is_in_list(g3, unique_list):
                    unique_list.append(g3)

for i, g in enumerate(unique_list):
    if not is_in_list(g, g3s):
        unique_list[i] = (g, "missed by Schwentner, et al.")
    else:
        unique_list[i] = (g, "constructed both ways")

write_graphs_to_file(unique_list, "arrow_merge_%d_to_%d"%(d1,d2))
