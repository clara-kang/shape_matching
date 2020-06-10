import bpy
import bmesh
import sys
import os


mat = bpy.data.materials.get("mat_single_pt")
if not mat:
    mat = bpy.data.materials.new("mat_single_pt")
    mat.diffuse_color = (1.0, 0.647, 0)

bpyscene = bpy.context.scene

def mark_v(loc, id):
    mesh = bpy.data.meshes.new('Basic_Sphere')
    basic_sphere = bpy.data.objects.new("v_" + str(id), mesh)
    bpyscene.objects.link(basic_sphere)
    bpyscene.objects.active = basic_sphere
    basic_sphere.select = True
    basic_sphere.location = loc
    basic_sphere.active_material = mat
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=1, diameter=0.02)
    bm.to_mesh(mesh)
    bm.free()

def mark_edge(edge, id):
    verts = edge.verts
    loc = ( verts[0].co + verts[1].co ) / 2
    mark_v(loc, id)

# Get the active mesh
me = bpy.context.object.data

patch_edges = {}
patch_ordr_edges = {}
patch_ordr_verts = {} # the verts are the end of edges
boundaries = []
bndry_patch_pairs = []
edge_assigned = {}
bndry_verts = []
bndry_verts_ids = {}
shape_curve_ids = []
# vid where a bndry starts [[start1, start2, start3], [start1, start2, start3, start4],...]
bndry_start_ids = {}

# Get a BMesh representation
bm = bmesh.new()   # create an empty BMesh
bm.from_mesh(me)   # fill it in from a Mesh

def insertToPatchEdges(idx, data):
    global patch_edges
    if idx in patch_edges:
        patch_edges[idx].add(data)
    else:
        patch_edges[idx] = {data}

# the boundaries for each patch, unordered
for e in bm.edges:
    if e.seam:
        # print("e: ", e)
        for adjFace in e.link_faces:
            insertToPatchEdges(adjFace.material_index, e)
            edge_assigned[e] = False

def getNextEdge(edges, last_e, v):
    link_es = v.link_edges
    for e in link_es:
        if e in edges and not (e == last_e):
            nv = e.verts[0]
            if nv == v:
                nv = e.verts[1]
            return (e, nv)

# order the edges for each patch
for pid in patch_edges:
    p_edges = patch_edges[pid]
    start_e = next(iter(p_edges))
    start_v = start_e.verts[0]
    ordered_edges = [start_e]
    ordered_verts = [start_v]
    next_e, next_v = getNextEdge(p_edges, start_e, start_v)
    while next_e != start_e:
        ordered_edges.append(next_e)
        ordered_verts.append(next_v)
        # print("next_v: ", next_v)
        next_e, next_v = getNextEdge(p_edges, next_e, next_v)
    # print("ordered_verts: ", ordered_verts)
    patch_ordr_edges[pid] = ordered_edges
    patch_ordr_verts[pid] = ordered_verts
# print("patch_ordr_edges: ", patch_ordr_edges)
# pid = 0
# for id in range(len(patch_ordr_edges[pid])):
#     mark_edge(patch_ordr_edges[pid][id], id)
# for id in range(len(patch_ordr_verts[pid])):
#     mark_v(patch_ordr_verts[pid][id].co, id)

# get boundry pairs

def get_adj_patches(faces):
    adj_patches = set()
    for face in faces:
        adj_patches.add(face.material_index)
    return adj_patches

# get a start v/e id for bndry
def getBndryStart(pid, start_eid):
    global patch_ordr_verts, patch_ordr_edges
    eid = start_eid
    start_edge = patch_ordr_edges[pid][eid]

    adj_patches = get_adj_patches(start_edge.link_faces)
    # print("adj_patches: ", adj_patches)
    while True:
        eid = ( eid + 1 ) % len(patch_ordr_edges[pid])
        next_edge = patch_ordr_edges[pid][eid]
        next_adj_patches = get_adj_patches(next_edge.link_faces)
        # print("next_adj_patches: ", next_adj_patches)
        if next_adj_patches != adj_patches:
            return ( eid - 1 + len(patch_ordr_edges[pid]) ) % len(patch_ordr_edges[pid])
        else:
            adj_patches = next_adj_patches

# get bndry vertices
uv_lay = bm.loops.layers.uv.active

def getVertUV(pid, vid):
    global patch_ordr_edges, patch_ordr_verts
    e = patch_ordr_edges[pid][vid]
    link_faces = e.link_faces
    for face in link_faces:
        # get the face inside the patch
        if face.material_index == pid:
            for loop in face.loops:
                # get uv of corresponding vert
                if loop.vert == patch_ordr_verts[pid][j]:
                    return loop[uv_lay].uv

pids = []

for pid in patch_ordr_edges:
    pids.append(pid)

# get bndry_start_ids
for pid in pids:
    start = getBndryStart(pid, 0)
    bndry_start_ids[pid] = [start]
    eid = ( start + 1 ) % (len(patch_ordr_edges[pid]))
    while True:
        next_start = getBndryStart(pid, eid)
        if next_start == start:
            break
        bndry_start_ids[pid].append(next_start)
        eid = ( next_start + 1 ) % (len(patch_ordr_edges[pid]))

# print("bndry_start_ids: ", bndry_start_ids)
# # test patch_ordr_edges
# for pid in pids:
#     for start in bndry_start_ids[pid]:
#         mark_v(patch_ordr_verts[pid][start].co, "1")

for pid in pids:
    for i in range(len(bndry_start_ids[pid])):
        start_vid = bndry_start_ids[pid][i]
        next_start_vid = bndry_start_ids[pid][(i + 1)%(len(bndry_start_ids[pid]))]
        this_bndry_verts = []
        j = start_vid
        while True:
            this_bndry_verts.append(getVertUV(pid, j))
            j = (j + 1) % len(patch_ordr_verts[pid])
            if j == next_start_vid:
                break
        bndry_verts.append(this_bndry_verts)

# get shape_curve_ids
shape_curve_ids = []
path_pairs = []

cnt = 0
for i in range(len(pids)):
    shape_curve_ids.append([])
    for j in range(len(bndry_start_ids[pid])):
        shape_curve_ids[-1].append(cnt)
        cnt += 1

def getConnPid(pid, eid):
    global patch_ordr_edges
    e = patch_ordr_edges[pid][eid]
    adj_patches = get_adj_patches(e.link_faces)
    for adj_pid in adj_patches:
        if adj_pid != pid:
            return adj_pid
    return -1

def getBndryId(pid, bndry_id):
    cnt = 0
    for id in pids:
        if id != pid:
            cnt += len(shape_curve_ids[pid])
        else:
            break
    return cnt + bndry_id

def whichBndry(pid, v):
    global patch_ordr_edges, bndry_start_ids

    vid = patch_ordr_verts[pid].index(v)
    # between which two
    for i in range(len(bndry_start_ids[pid])):
        start_id = bndry_start_ids[pid][i]
        next_start_id = bndry_start_ids[pid][(i + 1)%(len(bndry_start_ids[pid]))]
        print("start_id: ", start_id)
        print("next_start_id: ", next_start_id)
        print("vid: ", vid)
        if next_start_id > start_id:
            if vid > start_id and vid < next_start_id:
                same_dir = ( vid - start_id == 1 )
                return ( i, same_dir )
        else:
            if vid > start_id:
                same_dir = ( vid - start_id == 1 )
                return ( i, same_dir )
            elif vid < next_start_id:
                same_dir = ( vid + len(bndry_start_ids[pid]) - start_id == 1 )
                return ( i, same_dir )

print("bndry_start_ids: ", bndry_start_ids)

def path_used(bndry_id):
    global path_pairs
    for pair in path_pairs:
        if pair[0] == bndry_id or pair[1] == bndry_id:
            return True
    return False

# # get bndry pairs
cnt = 0
for pid in pids:
    for i in range(len(bndry_start_ids[pid])):
        bndry_id = getBndryId(pid, i)
        if path_used(bndry_id):
            continue
        # get the connecting patch
        second_vid = ( bndry_start_ids[pid][i] + 1 ) % (len(patch_ordr_verts[pid]))# also second vid
        # mark_v(patch_ordr_verts[pid][second_vid].co, "1")
        adj_pid = getConnPid(pid, second_vid)
        res = whichBndry(adj_pid, patch_ordr_verts[pid][second_vid])
        adj_bndry_id = getBndryId(adj_pid, res[0])
        pair_info = (bndry_id, adj_bndry_id, res[1])
        path_pairs.append(pair_info)

print("path_pairs: ", path_pairs)
print("pids: ", pids)
# mark_v(patch_ordr_verts[pid][start_bndry_id].co, pid)
for id in range(len(bndry_verts)):
    print("bndry lens: ", len(bndry_verts[id]))

bm.free()

# ================
import shape_info_lib
import dsplc_bndry_lib
import order_calc_lib
import write_lib

shape_info = shape_info_lib.ShapesInfo(shape_curve_ids, bndry_verts, path_pairs)
dsplc_bndry_calc = dsplc_bndry_lib.DsplcCalc(shape_info)
dsplc_bndry_calc.get_pair_bndry_dsplcmnt()
dsplc_bndry_calc.get_var_ids()
dsplc_bndry_calc.group_vids()
dsplc_bndry_calc.get_closed_groups()
dsplc_bndry_calc.group_vars()
dsplc_bndry_calc.buildGroupConnGraph()
dsplc_bndry_calc.getOrder()
dsplc_bndry_calc.calcVariables()
dsplc_bndry_calc.itrpRestOfCurves()

write_lib.write_pts_out(shape_info, "bndry_pts.txt")
write_lib.write_uvs_out(shape_info, "bndry_pts_uv.txt")
write_lib.write_info_out(shape_info, dsplc_bndry_calc, "info.txt", seg_len, shape_verts_num)
