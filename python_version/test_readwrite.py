from readwrite import *
import ThinSincereQuivers as tsq

K4 = tsq.ToricQuiver([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
cs = tsq.coneSystem(K4)

fname = "test_json_writing.json"
write_arrays_to_file(cs, fname)

c2s = read_arrays_from_file(fname)

cs = [x.transpose() for x in cs]
fps = [tsq.flowPolytope(K4, x).transpose() for x in tsq.referenceThetas(cs)]
write_arrays_to_file(fps, "flow_polytope.json")
