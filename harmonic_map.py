import bpy
import bmesh
import numpy as np
import math

# bpy.data.objects['Sphere'].select = True
# bpy.context.scene.objects.active = bpy.data.objects['Sphere']
# me = bpy.data.objects['Sphere'].data
me = bpy.context.object.data

# Get a BMesh representation
# bm = bmesh.new()   # create an empty BMesh
bm = bmesh.from_edit_mesh(me)   # fill it in from a Mesh

# gather verts on boundaries
verts_on_bndrys = set()

for e in bm.edges:
    if e.seam or e.is_boundary:
        verts_on_bndrys.add(e.verts[0])
        verts_on_bndrys.add(e.verts[1])

N = len(bm.verts) - len(verts_on_bndrys)
M = np.zeros([N, N])
rhs_us = np.zeros(N)
rhs_vs = np.zeros(N)

print("N: ", N)
# vertex to index map
vert_2_index = {}
index_2_vert = [None] * N

cnt = 0
for v in bm.verts:
    if not (v in verts_on_bndrys):
        index_2_vert[cnt] = v
        vert_2_index[v] = cnt
        cnt += 1

def getVnb(v, nb_edges):
    nb_vs = []
    for i in range(len(nb_edges)):
        e = nb_edges[i]
        for i in range(2):
            if e.verts[i] != v:
                nb_vs.append(e.verts[i])
    return nb_vs

uv_lay = bm.loops.layers.uv.active

def getVal(pid, v):
    nb_fs = v.link_faces
    for nb_f in nb_fs:
        if nb_f.material_index == pid:
            for loop in nb_f.loops:
                if loop.vert == v:
                    return loop[uv_lay].uv

def setUV(v, uv):
    nb_fs = v.link_faces
    for nb_f in nb_fs:
        for loop in nb_f.loops:
            if loop.vert == v:
                loop[uv_lay].uv[0] = uv[0]
                loop[uv_lay].uv[1] = uv[1]
def getPid(v):
    nb_fs = v.link_faces
    return nb_fs[0].material_index

def getWeight(e):
    weight = 0
    nb_fs = e.link_faces
    for nb_f in nb_fs:
        for loop in nb_f.loops:
            if loop.vert not in e.verts:
                angle = loop.calc_angle()
                weight += 1 / math.tan(angle)
    return weight

def getNormalizedWeights(nb_edges):
    weights = []
    for nb_id in range(len(nb_edges)):
        weight = getWeight(nb_edges[nb_id])
        weights.append(weight)
    weights = np.array(weights)
    weights = weights / np.sum(weights)
    return weights

# build M
for i in range(N):
    v = index_2_vert[i]
    nb_edges = v.link_edges
    nb_vs = getVnb(v, nb_edges)
    weights = getNormalizedWeights(nb_edges)
    # if i == 0:
    #     print("len(nb_vs)", len(nb_vs))
    pid = getPid(v)

    # self
    M[i][i] = 1
    for nb_id in range(len(nb_vs)):
        nb_v = nb_vs[nb_id]
        weight = weights[nb_id]
        print("weight: ", weight)
        if not (nb_v in verts_on_bndrys):
            nb_idx = vert_2_index[nb_v]
            # print("nb_idx: ", nb_idx)
            M[i][nb_idx] -= weight
        else:
            val = getVal(pid, nb_v)
            rhs_us[i] += weight * val[0]
            rhs_vs[i] += weight * val[1]

us = np.linalg.solve(M, rhs_us)
vs = np.linalg.solve(M, rhs_vs)

# set uv
for i in range(N):
    v = index_2_vert[i]
    uv = [us[i], vs[i]]
    # uv = [0, 0]
    setUV(v, uv)
