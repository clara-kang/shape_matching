import numpy as np

# trsfrm_funcs = [(False, np.array([1, 1])), (True, np.array([1, 1])), (True, np.array([-1, -1])), [False, np.array([-1, -1])]]
trsfrm_funcs = [(False, np.array([1, 1])), (True, np.array([-1, 1])), (True, np.array([1, -1])), [False, np.array([-1, -1])]]

def applyTrnsfrm(coord, func):
    trnsfrmd_coord = np.zeros(2)
    if func[0]:
        trnsfrmd_coord[0] = coord[1]
        trnsfrmd_coord[1] = coord[0]
    else:
        trnsfrmd_coord[0] = coord[0]
        trnsfrmd_coord[1] = coord[1]
    trnsfrmd_coord = np.multiply(trnsfrmd_coord, func[1])
    return trnsfrmd_coord
