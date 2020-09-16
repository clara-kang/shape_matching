import numpy as np
import triangle as tr
import harmonic_mapping
import scipy.optimize as opt
import math

class Singularity:
    def __init__(self, loc, vs, fgs, fid):
        self.loc = loc
        self.vs = vs
        self.bs = np.zeros(3)
        self.cs = np.zeros(3)

        self.bs[0] = vs[1][1] - vs[2][1]
        self.bs[1] = vs[2][1] - vs[0][1]
        self.bs[2] = vs[0][1] - vs[1][1]

        self.cs[0] = vs[2][0] - vs[1][0]
        self.cs[1] = vs[0][0] - vs[2][0]
        self.cs[2] = vs[1][0] - vs[0][0]

        self.twoA = np.dot(vs[:,0], self.bs)

        fs = fgs[:,0]
        gs = fgs[:,1]

        self.fid = fid

        self.fx = np.dot(fs, self.bs) / self.twoA
        self.fy = np.dot(fs, self.cs) / self.twoA
        self.gx = np.dot(gs, self.bs) / self.twoA
        self.gy = np.dot(gs, self.cs) / self.twoA

    def __repr__(self):
        return "fx: " + str(self.fx) + ", fy:" + str(self.fy) + ", gx:" + str(self.gx) + ", gy:" + str(self.gy)

    def fun(self, x):
        return (4*np.sin(x)*np.cos(x)**3-4*np.sin(x)**3*np.cos(x)) / (np.sin(x)**4+np.cos(x)**4-6*np.sin(x)**2*np.cos(x)**2) - (self.gx*np.cos(x)+self.gy*np.sin(x))/(self.fx*np.cos(x)+self.fy*np.sin(x))

    def find_sep_angles(self):
        guesses = np.linspace(0, 2 * np.pi, num=100, endpoint=False)
        sols = opt.fsolve(self.fun, guesses, xtol=1e-12)
        prc_sols = np.sort(np.unique(sols.round(decimals=8)))
        prc_sols_2 = []

        for x in prc_sols:
            if x < 0 or x > 2 * np.pi:
                continue
            left_angle = np.arctan2(np.sin(4*x), np.cos(4*x))
            right_angle = np.arctan2((self.gx*np.cos(x)+self.gy*np.sin(x)), (self.fx*np.cos(x)+self.fy*np.sin(x)))
            if round(left_angle, 4) == round(right_angle, 4):
                prc_sols_2.append(x)

        self.angles = prc_sols_2

class BoundaryNode:
    def __init__(self, loc, v_before):
        self.loc = loc
        self.v_before = v_before

class MegaNode:
    def __init__(self, loc, out_line_ids=[], list_id = 0):
        self.loc = loc
        self.out_line_ids = out_line_ids
        self.list_id = list_id

def quadrangulate(bndry_verts, shape_curve_ids, shape_id_test, bndry_reps):
    def segIntrsct(s1, e1, s2, e2):
        v1 = e1 - s1
        v2 = e2 - s2
        M = np.zeros([2,2])
        M[:,0] = v1
        M[:,1] = -v2
        rhs = s2 - s1
        ts = np.linalg.solve(M, rhs)
    #     print("ts: ", ts)
        if ts[0] >= 0 and ts[0] <= 1 and ts[1] >= 0 and ts[1] <= 1:
            return (True, s1 + v1 * ts[0])
        return (False, None)

    def getBaryCoords(v, vs):
        M = np.zeros([2,2])
        M[0] = vs[0] - vs[2]
        M[1] = vs[1] - vs[2]
        rhs = v - vs[2]
        bary_coords = np.linalg.solve(M.T, rhs)
        bary_coords = np.append(bary_coords, 1-bary_coords[0]-bary_coords[1])
        return bary_coords

    def inRange(coords):
        res = True
        for i in range(3):
            res = res and (coords[i] >= 0 and coords[i] <=1)
        return res

    def findFace(v):
        nonlocal faces, vertices

        for fid in range(len(faces)):
            vs = vertices[faces[fid]]
            bary_coords = getBaryCoords(v, vs)

            if inRange(bary_coords):
                return (fid, bary_coords)
        return None

    stream_lines = []

    def getClosestArm(rep, v_angle):
        rep_angle = np.arctan2(rep[1], rep[0])
        start_angle = rep_angle / 4
    #     v_angle = np.arctan2(v[1], v[0])
        v_angle = ( v_angle + 2 * np.pi ) % ( 2 * np.pi )
        min_diff = np.inf
        min_arm_angle = 0

        for i in range(4):
            arm_angle = ( start_angle + i * np.pi / 2 + 2 * np.pi ) % ( 2 * np.pi )
            diff = min( abs(arm_angle - v_angle), 2 * np.pi - abs(arm_angle - v_angle) )
            if diff < min_diff:
                min_diff = diff
                min_arm_angle = arm_angle

        return min_arm_angle

    def getClosestEdge(coords, fid):
        nonlocal faces, vertices
        smallest_dim = np.argmin(coords)
        if smallest_dim == 0:
            edge_vids = [1, 2]
        elif smallest_dim == 1:
            edge_vids = [0, 2]
        else:
            edge_vids = [0, 1]
        return faces[fid][edge_vids]

    def traceSeparatrix(src_sing, angle_id):
        nonlocal dt_reg, dt_init, stream_lines, vertices, faces

        def createBndryNode(v, fid, stream_line):
            edge = getClosestEdge(getBaryCoords(stream_line[-1], vertices[faces[fid]]), fid)
            edge_vs = vertices[edge]
            bndry_node_loc = segIntrsct(stream_line[-1], v, edge_vs[0], edge_vs[1])[1]

            # maybe into another triangle?
            if type(bndry_node_loc) is not np.ndarray:
                # get v on bndry:
                vid_on_bndry = min(faces[fid])
                bndry_node_loc = segIntrsct(stream_line[-1], v, vertices[vid_on_bndry], vertices[vid_on_bndry+1])[1]
                stream_line.append(bndry_node_loc)
                bndry_node = BoundaryNode(bndry_node_loc, vid_on_bndry)
                bndry_nodes.append(bndry_node)
                return

            stream_line.append(bndry_node_loc)
            bndry_node = BoundaryNode(bndry_node_loc, min(edge))
            bndry_nodes.append(bndry_node)

        stream_line = []
        stream_line.append(src_sing.loc)

        v = src_sing.loc + dt_init * np.array([np.cos(src_sing.angles[angle_id]), np.sin(src_sing.angles[angle_id])])
        res = findFace(v)
        fid = res[0]
        coords = res[1]
        stream_dir = np.array([np.cos(src_sing.angles[angle_id]), np.sin(src_sing.angles[angle_id])])
        arm_angle = src_sing.angles[angle_id]

        dt = dt_init
        while True:
            # stop if meeting another singularity
            met_sing = False
            for sing in singularities:
                if sing == src_sing:
                    continue
                if fid == sing.fid:
                    stream_line.append(v)
                    stream_line.append(sing.loc)
                    met_sing = True
                    break
            if met_sing:
                break


            face = faces[fid]
            stream_line.append(v)

            fgs = np.array(rep_vs[face])
            fg = coords[0] * fgs[0] + coords[1] * fgs[1] + coords[2] * fgs[2]
            fg_norm = np.linalg.norm(fg)

            arm_angle = getClosestArm(fg, arm_angle)
            arm_v = np.array([np.cos(arm_angle), np.sin(arm_angle)]) * fg_norm

            v = v + dt * arm_v

            test_coords = getBaryCoords(v, vertices[face])
            if inRange(test_coords):
                coords = test_coords
                continue
            else:
                dt = dt_reg
                res = findFace(v)
                if res == None:
                    createBndryNode(v, fid, stream_line)
                    break
                fid = res[0]
                coords = res[1]

        stream_line = np.array(stream_line)
        stream_lines.append(stream_line)

    def getCorners():
        nonlocal corners
        cnt = 0
        for i in range(len(shape_curve_ids[shape_id_test])):
            corners.append(cnt)
            cid = shape_curve_ids[shape_id_test][i]
            cnt += len(bndry_verts[cid])

    def streamLineIntrsctDumb(sid1, sid2):
        nonlocal stream_lines
        v_num_1 = len(stream_lines[sid1])
        v_num_2 = len(stream_lines[sid2])
        print("v_num_1: ", v_num_1)
        print("v_num_2: ", v_num_2)

        for i in range(v_num_1 - 1):
            for j in range(v_num_2 - 1):
                res = segIntrsct(stream_lines[sid1][i], stream_lines[sid1][i+1], stream_lines[sid2][j], stream_lines[sid2][j+1])
                if res[0]:
                    print(res[1])
                    end1 = i + 1
                    end2 = j + 1
                    intrsctn = res[1]

                    new_stream_line_1 = stream_lines[sid1][end1:v_num_1]
                    new_stream_line_2 = stream_lines[sid2][end2:v_num_2]
                    stream_lines[sid1] = stream_lines[sid1][:end1]
                    stream_lines[sid2] = stream_lines[sid2][:end2]

                    new_stream_line_1 = np.insert(new_stream_line_1, 0, intrsctn, axis=0)
                    new_stream_line_2 = np.insert(new_stream_line_2, 0, intrsctn, axis=0)
                    stream_lines[sid1] = np.append(stream_lines[sid1], [intrsctn], axis=0)
                    stream_lines[sid2] = np.append(stream_lines[sid2], [intrsctn], axis=0)

                    stream_lines.append(new_stream_line_1)
                    stream_lines.append(new_stream_line_2)
                    return

    def streamLineIntrsct(sid1, sid2):
        nonlocal stream_lines
        v_num_1 = len(stream_lines[sid1])
        v_num_2 = len(stream_lines[sid2])
        res = streamLineIntrsctRec(sid1, sid2, 0, v_num_1-1, 0, v_num_2-1)
        if res != None:
            (start1, end1, start2, end2, intrsctn) = res
            if start1==0 or start2==0:
                return
            print("intrsctn: ", intrsctn)
            display_intrsctns.append(intrsctn)

            new_stream_line_1 = stream_lines[sid1][end1:v_num_1]
            new_stream_line_2 = stream_lines[sid2][end2:v_num_2]
            stream_lines[sid1] = stream_lines[sid1][:end1]
            stream_lines[sid2] = stream_lines[sid2][:end2]

            new_stream_line_1 = np.insert(new_stream_line_1, 0, intrsctn, axis=0)
            new_stream_line_2 = np.insert(new_stream_line_2, 0, intrsctn, axis=0)
            stream_lines[sid1] = np.append(stream_lines[sid1], [intrsctn], axis=0)
            stream_lines[sid2] = np.append(stream_lines[sid2], [intrsctn], axis=0)

            stream_lines.append(new_stream_line_1)
            stream_lines.append(new_stream_line_2)

    def streamLineIntrsctRec(sid1, sid2, start1, end1, start2, end2):
        nonlocal stream_lines

        v_num_1 = len(stream_lines[sid1])
        v_num_2 = len(stream_lines[sid2])
        mid1 = int((start1 + end1) / 2)
        mid2 = int((start2 + end2) / 2)

        print("start1, ", start1, ", mid1, ", mid1, ", end1, ", end1)
        print("start2, ", start2, ", mid2, ", mid2, ", end2, ", end2)

        if end1 - start1 == 1 and end2 - start2 == 1:
            intersection = segIntrsct(stream_lines[sid1][start1], stream_lines[sid1][end1], stream_lines[sid2][start2], stream_lines[sid2][end2])[1]
            return (start1, end1, start2, end2, intersection)
        elif end1 - start1 == 1:
            # intersect first
            if segIntrsct(stream_lines[sid1][start1], stream_lines[sid1][end1], stream_lines[sid2][start2], stream_lines[sid2][mid2])[0]:
                return streamLineIntrsctRec(sid1, sid2, start1, end1, start2, mid2)
            # intersect second
            elif segIntrsct(stream_lines[sid1][start1], stream_lines[sid1][end1], stream_lines[sid2][mid2], stream_lines[sid2][end2])[0]:
                return streamLineIntrsctRec(sid1, sid2, start1, end1, mid2, end2)
            else:
                return None
        elif end2 - start2 == 1:
            # second intersect first
            if segIntrsct(stream_lines[sid1][start1], stream_lines[sid1][mid1], stream_lines[sid2][start2], stream_lines[sid2][end2])[0]:
                return streamLineIntrsctRec(sid1, sid2, start1, mid1, start2, end2)
            # second intersect second
            elif segIntrsct(stream_lines[sid1][mid1], stream_lines[sid1][end1], stream_lines[sid2][start2], stream_lines[sid2][end2])[0]:
                return streamLineIntrsctRec(sid1, sid2, mid1, end1, start2, end2)
            else:
                return None
        else:
            # first intersect first
            if segIntrsct(stream_lines[sid1][start1], stream_lines[sid1][mid1], stream_lines[sid2][start2], stream_lines[sid2][mid2])[0]:
                return streamLineIntrsctRec(sid1, sid2, start1, mid1, start2, mid2)
            # first intersect second
            if segIntrsct(stream_lines[sid1][start1], stream_lines[sid1][mid1], stream_lines[sid2][mid2], stream_lines[sid2][end2])[0]:
                return streamLineIntrsctRec(sid1, sid2, start1, mid1, mid2, end2)
            # second intersect first
            if segIntrsct(stream_lines[sid1][mid1], stream_lines[sid1][end1], stream_lines[sid2][start2], stream_lines[sid2][mid2])[0]:
                return streamLineIntrsctRec(sid1, sid2, mid1, end1, start2, mid2)
            # second intersect second
            if segIntrsct(stream_lines[sid1][mid1], stream_lines[sid1][end1], stream_lines[sid2][mid2], stream_lines[sid2][end2])[0]:
                return streamLineIntrsctRec(sid1, sid2, mid1, end1, mid2, end2)
            else:
                return None

    def getMegaNodes():
        nonlocal mega_nodes, all_stream_lines
        for lid in range(len(all_stream_lines)):
            stream_line = all_stream_lines[lid]
            loc1 = stream_line[0]
            loc2 = stream_line[-1]
            locs = [loc1, loc2]
            for i in range(2):
                key = (locs[i][0], locs[i][1])
                if key in mega_nodes:
                    mega_nodes[key].out_line_ids.append(lid)
                else:
                    m_node = MegaNode(locs[i], [lid])
                    mega_nodes[key] = m_node
        cur_id = 0
        for key in mega_nodes:
            mega_nodes_list.append(mega_nodes[key])
            mega_nodes[key].list_id = cur_id
            cur_id += 1

    # 0 for start node, -1 for end node
    def getLineNode(stream_line, which_node):
        nonlocal mega_nodes
        return mega_nodes[(stream_line[which_node][0], stream_line[which_node][1])]

    def buildMegaMap():
        nonlocal mega_nodes_list, all_stream_lines, mega_map

        mega_map = np.zeros([len(mega_nodes_list), len(mega_nodes_list)]) - 1

        for lid in range(len(all_stream_lines)):
            stream_line = all_stream_lines[lid]
            node_0 = getLineNode(stream_line, 0)
            node_1 = getLineNode(stream_line, -1)
            mega_map[node_0.list_id][node_1.list_id] = lid
            mega_map[node_1.list_id][node_0.list_id] = lid

    def getOtherEnd(stream_line, m_node):
        if np.array_equal(stream_line[0], m_node.loc):
            return getLineNode(stream_line, -1)
        else:
            return getLineNode(stream_line, 0)

    def orderOutLines():
        for m_node in mega_nodes_list:
            out_angles = [None] * len(m_node.out_line_ids)
            for i in range(len(m_node.out_line_ids)):
                lid = m_node.out_line_ids[i]
                othr_end = getOtherEnd(all_stream_lines[lid], m_node)
                vec = othr_end.loc - m_node.loc
                out_angles[i] = ( np.arctan2(vec[1], vec[0]), lid )
            out_angles.sort()
            for i in range(len(m_node.out_line_ids)):
                 m_node.out_line_ids[i] = out_angles[i][1]


    def getNextLineId(m_node, src_line_id):
        idx = m_node.out_line_ids.index(src_line_id)
        next_idx = (idx - 1 + len(m_node.out_line_ids)) % len(m_node.out_line_ids)
        return m_node.out_line_ids[next_idx]

    def whichEnd(stream_line, m_node):
        if np.array_equal(stream_line[0], m_node.loc):
            return 0
        else:
            return -1

    def getQuadLineIds():
        nonlocal quads_line_ids

        # skip boundary lines
        line_masks = [None] * len(stream_lines)
        for i in range(len(stream_lines)):
            line_masks[i] = np.array([False, False])

        for line_id in range(len(stream_lines)):
            ends = [0, -1]

            # check both directions
            for i in range(2):
                if line_masks[line_id][i] == True:
                    continue

                line_masks[line_id][i] = True

                end = ends[i]
                q_lids = [line_id]
                line_end_node = getLineNode(all_stream_lines[line_id], end)
                next_line_id = getNextLineId(line_end_node, line_id)

                while next_line_id != line_id:
                    q_lids.append(next_line_id)
                    line_end_node = getOtherEnd(all_stream_lines[next_line_id], line_end_node)

                    # mark is this line is one of streamline
                    if next_line_id < len(stream_lines):
                        if whichEnd(all_stream_lines[next_line_id], line_end_node) == 0:
                            line_masks[next_line_id][0] = True
                        else:
                            line_masks[next_line_id][1] = True

                    next_line_id = getNextLineId(line_end_node, next_line_id)
                quads_line_ids.append(q_lids)

    line_groups = []

    def getLineGroup():
        nonlocal line_groups, quads_line_ids
        for q_lids in quads_line_ids:
            s1 = {q_lids[0], q_lids[2]}
            s2 = {q_lids[1], q_lids[3]}

            s1_merged = False
            s2_merged = False

            for i in range(len(line_groups)):
                if len(line_groups[i].intersection(s1)) > 0:
                    line_groups[i] = line_groups[i].union(s1)
                    s1_merged = True

                if len(line_groups[i].intersection(s2)) > 0:
                    line_groups[i] = line_groups[i].union(s2)
                    s2_merged = True

            if not s1_merged:
                line_groups.append(s1)
            if not s2_merged:
                line_groups.append(s2)

    def streamLineLength(stream_line):
        line_len = 0
        for i in range(len(stream_line)-1):
            line_len += np.linalg.norm(stream_line[i+1] - stream_line[i])
        return line_len

    def getAllLinelength():
        nonlocal all_stream_lines, lines_lens, inc_len
        for i in range(len(all_stream_lines)):
            lines_lens[i] = streamLineLength(all_stream_lines[i])
        print("lines_lens: ", lines_lens)
        inc_len = min(lines_lens)
        if shape_id_test == 2:
            inc_len /= 10

    def getAllLineUVlengths():
        nonlocal all_stream_lines, lines_lens, line_groups
        for group_lids in line_groups:
            group_line_lens = np.zeros(len(group_lids))
            group_lids_list = list(group_lids)
            for i in range(len(group_lids_list)):
                lid = group_lids_list[i]
                group_line_lens[i] = lines_lens[lid]
            avg_line_len = np.average(group_line_lens)

            inc_num = math.floor(avg_line_len / inc_len)

            for i in range(len(group_lids_list)):
                lid = group_lids_list[i]
                lines_uv_lens[lid] = inc_num

    def getLinePts():
        nonlocal line_pts, all_stream_lines, lines_uv_lens, lines_lens
        for lid in range(len(all_stream_lines)):
            stream_line = all_stream_lines[lid]
            uv_len = lines_uv_lens[lid]
            dl = lines_lens[lid] / uv_len

            line_pts[lid] = []

            # do not put vert on endpoints
            if uv_len <= 1:
                continue

            acc_len = 0
            for vid in range(len(stream_line) - 1):
                seg_len = np.linalg.norm(stream_line[vid+1] - stream_line[vid])
                acc_len += seg_len
                if acc_len >= dl:
                    excess_len = acc_len - dl
                    excess_t = excess_len / seg_len
                    pt = excess_t * stream_line[vid] + (1 - excess_t) * stream_line[vid+1]
                    line_pts[lid].append(pt)
                    acc_len = excess_len
                    if len(line_pts[lid]) >= uv_len - 1:
                        break
            line_pts[lid] = np.array(line_pts[lid])

    def getCommonEndpoint(stream_line_1, stream_line_2):
        for i in [0, -1]:
            for j in [0, -1]:
                if np.array_equal(stream_line_1[i], stream_line_2[j]):
                    return (stream_line_1[i], i)

    def getPatchVerts():
        nonlocal all_stream_lines, quads_line_ids, line_pts
        for q_lids in quads_line_ids:
            patch_verts = []
            for i in range(len(q_lids)):
                lid = q_lids[i]
                next_lid = q_lids[(i+1)%len(q_lids)]
                (ep, s1_idx) = getCommonEndpoint(all_stream_lines[lid], all_stream_lines[next_lid])
                if s1_idx == -1:
                    for pt in line_pts[lid]:
                        patch_verts.append(pt)
                else:
                    for pt in list(reversed(line_pts[lid])):
                        patch_verts.append(pt)
                patch_verts.append(ep)
            patches_verts.append(patch_verts)

    def triangulatePatch(patch_id):
        segments = []
        for i in range(len(patches_verts[patch_id])-1):
            segments.append([i, i+1])
        segments.append([len(patches_verts[patch_id])-1, 0])

        A = dict(vertices=patches_verts[patch_id], segments=segments)
        B = tr.triangulate(A, 'pqa.02')
        return B

    def parameterizaPatch(trngl_info, patch_id):
        nonlocal line_pts, quads_line_ids
        q_lids = quads_line_ids[patch_id]
        uv_length = len(line_pts[q_lids[0]]) + 1
        uv_width = len(line_pts[q_lids[1]]) + 1

        uvs = []

        for i in range(1, uv_length+1):
            uvs.append(np.array([i, 0]))
        for i in range(1, uv_width + 1):
            uvs.append(np.array([uv_length, i]))
        for i in range(uv_length-1, -1, -1):
            uvs.append(np.array([i, uv_width]))
        for i in range(uv_width-1, -1, -1):
            uvs.append(np.array([0, i]))
        uvs = np.array(uvs)

        Us = harmonic_mapping.map(trngl_info['triangles'], trngl_info['vertices'], uvs[:,0])
        Vs = harmonic_mapping.map(trngl_info['triangles'], trngl_info['vertices'], uvs[:,1])
        UVs = np.stack((Us, Vs)).T
        return UVs


    def parameterizeAllPatches():
        file = open("objs/patches.obj", "w")

        def writeOut(trngl_info, uvs, patch_id, offset):

            file.write("o patch_"+str(patch_id)+"\n")
            for i in range(len(trngl_info['vertices'])):
                v = trngl_info['vertices'][i]
                file.write("v " + str(v[0]) + " " + str(v[1]) + " 0.0\n")
            file.write("vn 0.0 0.0 1.0\n")
            for i in range(len(uvs)):
                uv = uvs[i]
                file.write("vt " + str(uv[0]) + " " + str(uv[1]) + "\n")
            for face in trngl_info['triangles']:
                file.write("f ")
                for i in range(3):
                    file.write(str(face[i]+1+offset) + "/" + str(face[i]+1+offset) + "/" + str(patch_id+1) + " ")
                file.write("\n")

        offset = 0
        for patch_id in range(len(quads_line_ids)):
            trngl_info = triangulatePatch(patch_id)
            uvs = parameterizaPatch(trngl_info, patch_id)
            writeOut(trngl_info, uvs, patch_id, offset)
            offset += len(trngl_info['vertices'])
            print("len(trngl_info['vertices']): ", len(trngl_info['vertices']))

        file.close()

    # start ===================================================================
    # get boundary vertices of shape
    shape1_verts = []
    shape1_reps = []
    for cid in shape_curve_ids[shape_id_test]:
        for i in range(len(bndry_verts[cid])):
            shape1_verts.append(bndry_verts[cid][i])
            shape1_reps.append(bndry_reps[cid][i])

    # triangulation
    segments = []
    for i in range(len(shape1_verts)-1):
        segments.append([i, i+1])
    segments.append([len(shape1_verts)-1, 0])

    A = dict(vertices=shape1_verts, segments=segments)
    B = tr.triangulate(A, 'pqa.02')

    faces = B['triangles']
    vertices = B['vertices']
    bndry_n = len(shape1_verts)
    N = len(B['vertices'])

    # get representation angles
    Fs = harmonic_mapping.map(faces, vertices, np.array(shape1_reps)[:,0])
    Gs = harmonic_mapping.map(faces, vertices, np.array(shape1_reps)[:,1])
    rep_vs = np.stack((Fs, Gs)).T
    mem_as = [None]*len(rep_vs)

    # get member angles
    for i in range(len(rep_vs)):
        rep_angle = np.arctan2(rep_vs[i][1], rep_vs[i][0])
        angle = rep_angle / 4
        mem_as[i] = angle

    # find singularities
    singularities = []

    for fid in range(len(faces)):
        face = faces[fid]
        vs = vertices[face]
        fgs = rep_vs[face]
        M = np.array([[fgs[0][0] - fgs[2][0], fgs[1][0] - fgs[2][0]], [fgs[0][1] - fgs[2][1], fgs[1][1] - fgs[2][1]]])

        if np.linalg.matrix_rank(M) < 2:
            continue
        rhs = np.array([-fgs[2][0], -fgs[2][1]])
        sing_coords = np.linalg.solve(M, rhs)
        if 0 <= sing_coords[0] <= 1 and 0 <= sing_coords[1] <= 1 and 0 <= sing_coords[0] + sing_coords[1] <= 1:
            sing_pos = sing_coords[0] * vs[0] + sing_coords[1] * vs[1] + (1 - sing_coords[0] - sing_coords[1]) * vs[2]
            sing = Singularity(sing_pos, vs, fgs, fid)
            sing.find_sep_angles()
            singularities.append(sing)

    # streamlines
    dt_reg = 0.01
    dt_init = 0.01
    bndry_nodes = []
    stream_lines = []

    for sing in singularities:
        for a_id in range(len(sing.angles)):
            traceSeparatrix(sing, a_id)

    corners = []
    getCorners()
    corner_vs = vertices[corners]

    bndry_nodes.sort(key=lambda x: x.v_before)
    bndry_stream_lines = []

    for i in range(len(bndry_nodes)):
        stream_line_1 = [bndry_nodes[i].loc]
        stream_line_2 = []
        start_vid = bndry_nodes[i].v_before + 1
        stop_vid = bndry_nodes[(i+1)%len(bndry_nodes)].v_before+1
        stop_mid = False

        if stop_vid < start_vid:
            stop_vid += len(shape1_verts)

        # check if stream line contain corners
        for c in corners:
            if bndry_nodes[(i+1)%len(bndry_nodes)].v_before+1 < start_vid:
                if ( start_vid < c + len(shape1_verts) and c + len(shape1_verts) < stop_vid ):
                    stop_mid = True
                    mid = c + len(shape1_verts)
            if start_vid < c and c < stop_vid:
                stop_mid = True
                mid = c

        if stop_mid:
            for j in range(start_vid, mid+1):
                stream_line_1.append(shape1_verts[j % len(shape1_verts)])
            stream_line_1 = np.array(stream_line_1)
            bndry_stream_lines.append(stream_line_1)

            for j in range(mid, stop_vid):
                stream_line_2.append(shape1_verts[j % len(shape1_verts)])

            stream_line_2.append(bndry_nodes[(i+1)%len(bndry_nodes)].loc)
            stream_line_2 = np.array(stream_line_2)
            bndry_stream_lines.append(stream_line_2)

        else:
            for j in range(start_vid, stop_vid):
                stream_line_1.append(shape1_verts[j % len(shape1_verts)])
            stream_line_1.append(bndry_nodes[(i+1)%len(bndry_nodes)].loc)
            stream_line_1 = np.array(stream_line_1)
            bndry_stream_lines.append(stream_line_1)

    for i in range(5):
        for j in range(5, 10):
            print("i: ", i)
            print("j: ", j)
            streamLineIntrsctDumb(i, j)

    all_stream_lines = stream_lines + bndry_stream_lines

    mega_nodes = {}
    mega_nodes_list = []
    mega_map = None
    getMegaNodes()
    buildMegaMap()
    orderOutLines()

    quads_line_ids = []
    getQuadLineIds()

    line_groups = []
    getLineGroup()

    lines_uv_lens = [0] * len(all_stream_lines)
    lines_lens = [0] * len(all_stream_lines)
    inc_len = 0
    getAllLinelength()
    getAllLineUVlengths()

    line_pts = [None] * len(all_stream_lines)
    getLinePts()

    patches_verts = []
    getPatchVerts()
    parameterizeAllPatches()

    if shape_id_test == 2:
        bndry_stream_lines = [None] * 4
        print("shape_curve_ids[shape_id_test]: ", shape_curve_ids[shape_id_test])
        bndry_stream_lines[0] = bndry_verts[12]
        bndry_stream_lines[0] = np.append(bndry_stream_lines[0], [bndry_verts[13][0]], axis = 0)
        bndry_stream_lines[1] = np.append(bndry_verts[13], bndry_verts[14], axis = 0)
        bndry_stream_lines[1] = np.append(bndry_stream_lines[1], [bndry_verts[15][0]], axis = 0)
        bndry_stream_lines[2] = bndry_verts[15]
        bndry_stream_lines[2] = np.append(bndry_stream_lines[2], [bndry_verts[16][0]], axis = 0)
        bndry_stream_lines[3] = bndry_verts[16]
        bndry_stream_lines[3] = np.append(bndry_stream_lines[3], [bndry_verts[12][0]], axis = 0)
