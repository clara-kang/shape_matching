import bpy
import bmesh
import numpy as np

bpy.data.objects['Sphere'].select = True
bpy.context.scene.objects.active = bpy.data.objects['Sphere']
me = bpy.data.objects['Sphere'].data

# Get a BMesh representation
# bm = bmesh.new()   # create an empty BMesh
bm = bmesh.from_edit_mesh(me)   # fill it in from a Mesh

# gather verts on boundaries
verts_on_bndrys = set()

for e in bm.edges:
    if e.seam:
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

def getVnb(v):
    global bm
    nb_edges = v.link_edges
    nb_vs = set()
    for e in nb_edges:
        for i in range(2):
            if e.verts[i] != v:
                nb_vs.add(e.verts[i])
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

# build M
for i in range(N):
    v = index_2_vert[i]
    nb_vs = getVnb(v)
    if i == 0:
        print("len(nb_vs)", len(nb_vs))

    weight = 1 / len(nb_vs)
    pid = getPid(v)

    # self
    M[i][i] = 1
    for nb_v in nb_vs:
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
