import bpy
import bmesh

import sys
import os

# blend_dir = os.path.basename(bpy.data.filepath)
# print("blend_dir: ",  blend_dir)
# if blend_dir not in sys.path:
#    sys.path.append(blend_dir)
sys.path.append(r'C:\Users\Clara\iCloudDrive\circle_square_map\shape_matching')

import trnsfrm_func_lib
import shape_info_lib
import dsplc_bndry_lib
import order_calc_lib
import write_lib

import imp
imp.reload(shape_info_lib)
imp.reload(dsplc_bndry_lib)
imp.reload(order_calc_lib)
imp.reload(write_lib)
imp.reload(trnsfrm_func_lib)

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


# bpy.data.objects['Tee_S'].select = True
# bpy.context.scene.objects.active = bpy.data.objects['Tee_S']
# bpy.ops.object.mode_set(mode = 'OBJECT')

me = bpy.context.object.data

patch_edges = {}
patch_ordr_edges = {}
patch_ordr_verts = {} # the verts are the end of edges
boundaries = []
bndry_patch_pairs = []
edge_assigned = {}
bndry_verts = []
bndry_verts_ids = []
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
    if e.seam or e.is_boundary:
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

size_multiplier = 10
def getVertUV(pid, vid):
    global patch_ordr_edges, patch_ordr_verts, uv_lay
    e = patch_ordr_edges[pid][vid]
    link_faces = e.link_faces
    for face in link_faces:
        # get the face inside the patch
        if face.material_index == pid:
            for loop in face.loops:
                # get uv of corresponding vert
                if loop.vert == patch_ordr_verts[pid][vid]:
                    return [loop[uv_lay].uv[0] * size_multiplier, loop[uv_lay].uv[1] * size_multiplier]

def setVertUV(pid, vid, uv):
    global patch_ordr_edges, patch_ordr_verts
    v = patch_ordr_verts[pid][vid]
    link_faces = v.link_faces
    for face in link_faces:
        # get the face inside the patch
        if face.material_index == pid:
            for loop in face.loops:
                # get uv of corresponding vert
                if loop.vert == patch_ordr_verts[pid][vid]:
                    loop[uv_lay].uv[0] = uv[0]
                    loop[uv_lay].uv[1] = uv[1]

pids = []

for pid in patch_ordr_edges:
    pids.append(pid)

# for pid in pids:
#     print("len of ", pid, ": ", len(patch_ordr_verts[pid]))
# print("pids: ", pids)

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

# get bndry_verts
for pid in pids:
    for i in range(len(bndry_start_ids[pid])):
        start_vid = bndry_start_ids[pid][i]
        next_start_vid = bndry_start_ids[pid][(i + 1)%(len(bndry_start_ids[pid]))]
        this_bndry_verts = []
        this_bndry_vert_ids = []
        j = start_vid
        while True:
            this_bndry_verts.append(getVertUV(pid, j))
            this_bndry_vert_ids.append(j)
            j = (j + 1) % len(patch_ordr_verts[pid])
            if j == next_start_vid:
                break

        # make it straight
        # end_uv = getVertUV(pid, next_start_vid)
        # for i in range(len(this_bndry_verts)):
        #     t = i / len(this_bndry_verts)
        #     for dim in range(2):
        #         this_bndry_verts[i][dim] = (1 - t) * this_bndry_verts[0][dim] + t * end_uv[dim]
        bndry_verts.append(this_bndry_verts)
        bndry_verts_ids.append(this_bndry_vert_ids)

# for id in range(len(bndry_verts)):
#     print(id, " ", "bndry lens: ", len(bndry_verts[id]))

for id in range(len(bndry_verts_ids)):
    print(id, ": ", bndry_verts_ids[id])

# print("bndry_verts: ", bndry_verts)
# get shape_curve_ids
shape_curve_ids = []
path_pairs = []

# print("bndry_start_ids: ", bndry_start_ids)

cnt = 0
for pid in pids:
    shape_curve_ids.append([])
    for j in range(len(bndry_start_ids[pid])):
        shape_curve_ids[-1].append(cnt)
        cnt += 1

for pid in range(len(shape_curve_ids)):
    print( pid, ": ", shape_curve_ids[pid])
# print("shape_curve_ids: ", shape_curve_ids)

def getConnPid(pid, eid):
    global patch_ordr_edges
    e = patch_ordr_edges[pid][eid]
    adj_patches = get_adj_patches(e.link_faces)
    for adj_pid in adj_patches:
        if adj_pid != pid:
            return adj_pid
    return -1

def getBndryId(pid, which_bndry):
    cnt = 0
    for id in pids:
        if id != pid:
            cnt += len(shape_curve_ids[id])
        else:
            break
    return cnt + which_bndry

# def whichBndry(pid, v):
#     global shape_curve_ids, patch_ordr_verts
#     shape_bndry_ids = shape_curve_ids[pid]
#     vid = patch_ordr_verts[pid].index(v)
#     for bndry_id in shape_bndry_ids:
#         # print("bndry_id: ", bndry_id)
#         # print("bndry_verts_ids[bndry_id]: ", bndry_verts_ids[bndry_id])
#         if vid in bndry_verts_ids[bndry_id]:
#             # print("bndry_id: ", bndry_id)
#             same_dir = ( bndry_verts_ids[bndry_id].index(vid) == 2 )
#             return ( bndry_id, same_dir )

def whichBndry(pid, v):
    global patch_ordr_edges, bndry_start_ids

    vid = patch_ordr_verts[pid].index(v)
    # print("bndry_start_ids[pid]: ", bndry_start_ids[pid])
    # print("vid: ", vid)
    # between which two
    for i in range(len(bndry_start_ids[pid])):
        start_id = bndry_start_ids[pid][i]
        next_start_id = bndry_start_ids[pid][(i + 1)%(len(bndry_start_ids[pid]))]
        if next_start_id > start_id:
            if vid > start_id and vid < next_start_id:
                same_dir = ( vid - start_id == 1 )
                # print("i: ", i)
                return ( i, same_dir )
        else:
            if vid > start_id:
                same_dir = ( vid - start_id == 1 )
                # print("i: ", i)
                return ( i, same_dir )
            elif vid < next_start_id:
                same_dir = ( vid + len(bndry_start_ids[pid]) - start_id == 1 )
                # print("i: ", i)
                return ( i, same_dir )

# print("bndry_start_ids: ", bndry_start_ids)

def path_used(bndry_id):
    global path_pairs
    for pair in path_pairs:
        if pair[0] == bndry_id or pair[1] == bndry_id:
            return True
    return False

skip_times = 0
# # get bndry pairs
for pid in pids:
    for i in range(len(bndry_start_ids[pid])):
        bndry_id = getBndryId(pid, i)
        if path_used(bndry_id):
            continue
        # get the connecting patch
        second_vid = ( bndry_start_ids[pid][i] + 1 ) % (len(patch_ordr_verts[pid]))# also second vid
        # mark_v(patch_ordr_verts[pid][second_vid].co, "1")
        adj_pid = getConnPid(pid, second_vid)
        if adj_pid < 0:
            skip_times += 1
            continue
        print("adj_pid: ", adj_pid)
        res = whichBndry(adj_pid, patch_ordr_verts[pid][second_vid])
        adj_bndry_id = getBndryId(adj_pid, res[0])
        # pair_info = (bndry_id, res[0], res[1])
        pair_info = (bndry_id, adj_bndry_id, res[1])
        path_pairs.append(pair_info)

print("skip_times: ", skip_times)
print("path_pairs: ", path_pairs)

# for pid in [0, 9]:
#     for i in range(len(bndry_start_ids[pid])):
#         start_v = bndry_start_ids[pid][i]
#         next_v =  bndry_start_ids[pid][(i + 1) % (len(bndry_start_ids[0]))]
#         j = start_v
#         cnt = 0
#         while True:
#             mark_v(patch_ordr_verts[pid][j].co, str(pid)+"_"+str(i)+"_"+str(cnt))
#             cnt += 1
#             j = ( j + 1 ) % len(patch_ordr_verts[pid])
#             if j == next_v:
#                 break



# print("shape_curve_ids: ", shape_curve_ids)
bpy.ops.wm.save_as_mainfile(filepath=bpy.data.filepath)
# ================


shape_info = shape_info_lib.ShapesInfo(shape_curve_ids, bndry_verts, path_pairs)
write_lib.write_pts_out(shape_info, "bndry_pts.txt")
# print("shape_info.bndry_verts_flat: ", shape_info.bndry_verts_flat)
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

#=========
cnt = 0
# for pid in pids:
#     for i in range(len(patch_ordr_verts[pid])):
#         # print("before: ", shape_info.bndry_verts_flat[cnt])
#         # print("after: ", shape_info.bndry_vert_uv[cnt])
#         setVertUV(pid, i, shape_info.bndry_vert_uv[cnt])
#         cnt += 1

for pid in pids:
    start_vid = bndry_start_ids[pid][0]
    j = start_vid
    while True:
        setVertUV(pid, j, shape_info.bndry_vert_uv[cnt])
        cnt += 1
        j = (j + 1) % len(patch_ordr_verts[pid])
        if j == start_vid:
            break
write_lib.write_pts_out(shape_info, "bndry_pts.txt")
write_lib.write_uvs_out(shape_info, "bndry_pts_uv.txt")

bm.to_mesh(me)
bm.free()

bpy.ops.wm.save_as_mainfile(filepath=bpy.data.filepath)
print("done")
