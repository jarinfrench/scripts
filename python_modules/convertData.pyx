cimport numpy as np
import numpy as np
cpdef dict convertData(df, indices):
    cdef dict res = {i: [0,0,0] for i in indices}
    for i in df.itertuples(index = False):
        res[i.timestep, i.ring_num][0] += 1
        if i.type == 2:
            res[i.timestep, i.ring_num][1] += 1
        elif i.type == 3:
            res[i.timestep, i.ring_num][2] += 1
    return res
