import json
import numpy as np

# works for output of coneSyste(Q)
def write_arrays_to_file(list_of_arrays, filename="output.json"):
    data = [x.tolist() for x in list_of_arrays]
    with open(filename, "w") as f:
         json.dump(data, f, ensure_ascii=False, indent=4)

def read_arrays_from_file(filename):
    with open(filename, "r") as f:
        data = [np.array(x) for x in json.load(f)]
    return data
