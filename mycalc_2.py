from svg.path import parse_path
import xml.etree.ElementTree as ET
import numpy as np

ns = "{http://www.w3.org/2000/svg}"
tree = ET.parse("svgs/doll_dress.svg")
root = tree.getroot()

g = root.findall(ns + "g")[0]
paths = g.findall(ns  + "path")
path_objs = []

for p in paths:
    path_obj = parse_path(p.attrib["d"])
    path_objs.append(path_obj)


html_str_1 = '''
<!DOCTYPE html>
<html>
<body>

<svg height="800" width="1000">
<g transform="scale(2.5)">
'''

html_path_str = ""

# for p in paths:
#     html_path_str += "<path d=\"" + p.attrib["d"] + "\" />\n"

# [cps for curve1, cps for curve2, ...]
paths_cps = []
# [length of curve1, length of curve2, ...]
paths_lengths = []
# [[curve1_id, curve2_id, ..], [curve1_id, curve2_id, ..], ...]
shape_curve_ids = []

c_id = 0
for p_obj in path_objs:
    path_cps = []
    path_lens = []
    curve_ids = []
    for curve in p_obj:
        p_name = type(curve).__name__
        d_str = ""
        if p_name == "CubicBezier":
            cpts = np.zeros([4, 2])
            cpts[0][0] = curve.start.real
            cpts[0][1] = -curve.start.imag
            cpts[1][0] = curve.control1.real
            cpts[1][1] = -curve.control1.imag
            cpts[2][0] = curve.control2.real
            cpts[2][1] = -curve.control2.imag
            cpts[3][0] = curve.end.real
            cpts[3][1] = -curve.end.imag

            d_str = d_str + "M " + str(curve.start.real) + " " + str(curve.start.imag) + " "
            d_str = d_str + "C " + str(curve.control1.real) + " " + str(curve.control1.imag)
            d_str = d_str + ", " + str(curve.control2.real) + " " + str(curve.control2.imag)
            d_str = d_str + ", " + str(curve.end.real) + " " + str(curve.end.imag)

            path_cps.append(cpts)
            curve_ids.append(c_id)
        elif p_name == "Line" or p_name=="Close":
            if curve.length() == 0:
                continue
            cpts = np.zeros([2, 2])
            cpts[0][0] = curve.start.real
            cpts[0][1] = -curve.start.imag
            cpts[1][0] = curve.end.real
            cpts[1][1] = -curve.end.imag

            d_str = d_str + "M " + str(curve.start.real) + " " + str(curve.start.imag) + " "
            d_str = d_str + "L " + str(curve.end.real) + " " + str(curve.end.imag)

            path_cps.append(cpts)
            curve_ids.append(c_id)
        else:
            continue

        path_lens.append(curve.length())

        html_path_str += "<path id = \"" + str(c_id) + "\" d=\"" + d_str + "\" stroke=\"black\" fill=\"transparent\" />\n"
        c_id +=1

    shape_curve_ids.append(curve_ids)
    paths_cps.extend(path_cps)
    paths_lengths.extend(path_lens)

# label the paths
html_text_str = "<text font-size=\"7\">"

for i in range(len(paths_cps)):
    html_text_str = html_text_str + "<textPath href=\"#" + str(i) + "\" startOffset=\"40%\">" + str(i) + "</textPath>\n"
html_text_str += "</text>"

# end of html file
html_str_2 = '''
</g>
</svg>

</body>
</html>
'''
html_file = open("svg_display.html", "w")
html_file.write(html_str_1)
html_file.write(html_path_str)
html_file.write(html_text_str)
html_file.write(html_str_2)
html_file.close()

shape_ccw = [None] * len(shape_curve_ids)
for sid in range(len(shape_curve_ids)):
    cids = shape_curve_ids[sid]
    v1 = paths_cps[cids[0]][-1] - paths_cps[cids[0]][0]
    wind_num = 0

    def getAngle(v1, v2):
        cross_p = np.cross(v1, v2)
        cos_theta = np.dot(v1, v2);
        if cross_p > 0:
            return cos_theta
        else:
            return -cos_theta

    for cid in range(1, len(cids)):
        v2 = paths_cps[cids[cid]][-1] - paths_cps[cids[cid]][0]
        wind_num += getAngle(v1, v2)
        v1 = v2;

    v2 = paths_cps[cids[0]][-1] - paths_cps[cids[0]][0]
    wind_num += getAngle(v1, v2)

    if wind_num > 0:
        shape_ccw[sid] = True
    else:
        shape_ccw[sid] = False
print("shape_ccw: ", shape_ccw)

# for bodice
path_pairs_no_orientation = [
    (8, 25),
    (22, 5),
    (13, 28),
    (14, 19),
    (30, 11),
    (31, 2),
    (23, 18),
    (1, 6),
    (27, 20),
    (3, 10),
    (12, 15),
    (32, 29)
]

def getShapeId(cid):
    global shape_curve_ids
    for sid in range(len(shape_curve_ids)):
        if cid in shape_curve_ids[sid]:
            return sid

path_pairs = []
for ppair in path_pairs_no_orientation:
    sid1 = getShapeId(ppair[0])
    sid2 = getShapeId(ppair[1])
    ccw1 = shape_ccw[sid1]
    ccw2 = shape_ccw[sid2]
    path_pairs.append((ppair[0], ppair[1], ccw1 != ccw2))

import math

path_lens = np.array(path_lens)
print("path_lens: ", path_lens)
# at least have this many pts on one curve
min_curve_pts_num = 6
# desired segment length
seg_len = np.amin(path_lens) / min_curve_pts_num

# get number of vertices for each curve
curve_vert_nums = [0] * len(paths_cps)

for pair in path_pairs:
    # use the maximum curve length of the pair of curves
    c1_len = paths_lengths[pair[0]]
    c2_len = paths_lengths[pair[1]]
    max_len = max(c1_len, c2_len)
    # ensure that each segment is less than seg_len
    c_vert_num = math.ceil( max_len / seg_len)
    # set vertices num for both curves
    curve_vert_nums[pair[0]] = c_vert_num
    curve_vert_nums[pair[1]] = c_vert_num

for cid in range(len(paths_cps)):
    # skip if already set
    if curve_vert_nums[cid] == 0:
        c_len = paths_lengths[cid]
        curve_vert_nums[cid] = math.ceil( c_len / seg_len)

bndry_verts = [None] * len(paths_cps)
err_thresh = 0.01
max_it = 20

def get_pair(c_id):
    global path_pairs
    for pair in path_pairs:
        if c_id == pair[0]:
            return pair[1]
        elif c_id == pair[1]:
            return pair[0]
    return -1

def getPtOnCurve(cps, t):
    return (1-t)**3 * cps[0] + 3*t*(1-t)**2 * cps[1] + 3*(1-t)*t**2 * cps[2] + t**3 * cps[3]

def getEquitDistantPtOnCurve(cps, pt_num, curve_len):
    N = 500
    prev_pt = cps[0]
    seg_lens = np.zeros(N)
    for i in range(1, N):
        t = i / N
        pt = getPtOnCurve(cps, t)
        seg_lens[i] = np.linalg.norm(pt - prev_pt)
        prev_pt = pt
    avg_seg_len = curve_len / pt_num
    bndry_pts = np.zeros([pt_num, 2])
    bndry_pts[0] = cps[0]

    last_seg = 0
    last_seg_part_len = 0
    ts = np.zeros(pt_num)

    for pt_id in range(1, pt_num):
        # find the segment up to which has avg_seg_len
        acc_len = last_seg_part_len
        for i in range(last_seg, N):
            if acc_len + seg_lens[i] < avg_seg_len:
                acc_len += seg_lens[i]
            else:
                # find the portion in segment
                last_seg_part_len = acc_len + seg_lens[i] - avg_seg_len
                last_seg = i + 1
                t = i / N - ( last_seg_part_len / seg_lens[i] ) / N
                ts[pt_id] = t
                bndry_pts[pt_id] = getPtOnCurve(cps, t)
#                 print("seg len = ", np.linalg.norm(bndry_pts[pt_id] - bndry_pts[pt_id-1]))
                break
    return bndry_pts

for curve_id in range(len(paths_cps)):
    curve_cps = paths_cps[curve_id]
    # number of verts on curve already calculated and stored in curve_vert_nums
    bndry_vert_num = curve_vert_nums[curve_id]
    # store curve vertices here
    curve_verts = np.zeros([bndry_vert_num, 2])

    # bezier curve
    if (len(curve_cps) == 4):
        curve_verts = getEquitDistantPtOnCurve(curve_cps, bndry_vert_num, paths_lengths[curve_id])
    # line
    elif (len(curve_cps) == 2):
        avg_seg_len = paths_lengths[curve_id] / bndry_vert_num
        d_t = 1 / bndry_vert_num
        t = 0
        dir_vec = (curve_cps[1] - curve_cps[0])
        for v_id in range(bndry_vert_num):
            pt_loc = curve_cps[0] + t * dir_vec
            curve_verts[v_id] = pt_loc
            t += d_t
    bndry_verts[curve_id] = curve_verts

# resize the shape
size_multiplier = 1 / 16
for cid in range(len(bndry_verts)):
    for vid in range(len(bndry_verts[cid])):
        bndry_verts[cid][vid] *= size_multiplier

# how many vertices in each shape
shape_verts_num = [0] * len(path_objs)

for i in range(len(shape_curve_ids)):
    for curve_id in shape_curve_ids[i]:
        shape_verts_num[i] += len(bndry_verts[curve_id])


import get_polygons_2
import write_lib1
import subprocess
import numpy as np

polygon = get_polygons_2.Polygons(bndry_verts, shape_curve_ids, path_pairs)
polygon.getPolyVerts()
polygon.getPolyVecs()
polygon.getAngles()
polygon.getNormCenterSegDist()
polygon.getPairAngleDiff()
polygon.writeOut()
polygon.writeVertsOut("poly_verts.txt")
write_lib1.write_pts_out(polygon, "bndry_pts.txt")

subprocess.call(['C:\\gurobi902\\win64\\python37\\bin\\python.exe', 'opt_angles_5.py'])
polygon.readInRotation()
polygon.rotateShapes()
polygon.writeVertsOut("poly_verts_rot.txt")
polygon.getPolyVecs()
polygon.getAngles()
polygon.getPairAngleDiff()

polygon.writeOut2()
subprocess.call(['C:\\gurobi902\\win64\\python37\\bin\\python.exe', 'opt_verts_4.py'])

polygon.readInVerts()
polygon.intrpCurve()

bndry_c_verts_start = polygon.get_bndry_c_verts_start()

write_lib1.write_uvs_out(polygon, "bndry_pts_uv.txt")
write_lib1.write_info_out(polygon, "info.txt", seg_len, shape_verts_num, bndry_c_verts_start)
