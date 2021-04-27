import sys
sys.path.append("..")
import ToricQuiver as qc
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
    return gc.unoriented_unique_up_to_isomorphism(graphs)

d = int(input("Enter value for d -> "))

gs = read_step_file("../outputs/d=%d/step5"%d)
print("gs = ", gs)

for g in gs:
    max_unstables = qc.all_maximal_unstable(g)
    print(g, max_unstables)
    n = [len(m) for m in max_unstables] + [0]
    print("the quiver \n", g, "\nis ", g.shape[1] - max(n), " neighborly")

