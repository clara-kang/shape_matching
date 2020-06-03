from svg.path import parse_path
import xml.etree.ElementTree as ET
import numpy as np

ns = "{http://www.w3.org/2000/svg}"
tree = ET.parse("../bear_detached1.svg")
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

# for t-shirt
# path_pairs = [(3, 11, False), (7, 15, True), (19, 27, False), (23, 31, True), (8, 28, True), (10, 30, True), (24, 12, True), (26, 14, True), (0, 18, False), (2, 16, False), (6, 20, False), (4, 22, False)]
# for cube
# path_pairs = [(0, 1, False), (3, 4, False), (5, 2, False), (6, 13, False), (7, 10, False), (8, 9, False), (11, 12, False)]
# for cube detached
# path_pairs = [(0, 18, False), (1, 4, False), (2, 8, False), (3, 20, False), (5, 17, False), (6, 13, False), (7, 9, False), \
#     (10, 12, False), (11, 21, False), (14, 16, False), (15, 22, False), (19, 23, False)]
# for rotation test
# path_pairs = [(0, 4, False)]
# for bear
# path_pairs = [(0, 8, True), (1, 9, True), (4, 12, True), (5, 13, True), (14, 15, False), (10, 11, False), (2, 3, False), (6, 7, False)]
# for bear detached
# path_pairs = [(6, 16, True), (7, 17, True), (1, 11, True), (2, 12, True), (10, 15, True), \
#               (14, 19, True), (13, 18, True), (5, 0, True), (9, 4, True), (8, 3, True)]
# for bear detached 1
path_pairs = [(16, 20, False), (18, 24, False), (22, 27, False), (13, 8, False), (10, 4, False), \
    (6, 2, False), (0, 14, True), (5, 19, True), (9, 23, True), (12, 26, True), \
    (30, 34, False), (32, 38, False), (36, 41, False), (55, 50, False), (52, 46, False), (48, 44, False), \
    (42, 28, True), (47, 33, True), (51, 37, True), (54, 40, True), \
    (1, 29, True), (3, 31, True), (7, 35, True), (11, 39, True), (25, 53, True), (21, 49, True), (17, 45, True), \
    (15, 43, True)]

import math

path_lens = np.array(path_lens)
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


# how many vertices in each shape
shape_verts_num = [0] * len(path_objs)

for i in range(len(shape_curve_ids)):
    for curve_id in shape_curve_ids[i]:
        shape_verts_num[i] += len(bndry_verts[curve_id])


# pts_file = open("bndry_pts.txt", "w")
#
# for curve_id in range(len(paths_cps)):
#     for v in bndry_verts[curve_id]:
#         pts_file.write(str(v[0]) + " " + str(v[1]) + "\n")
#
# pts_file.close()


import shape_info_lib
import dsplc_bndry_lib
import order_calc_lib
import write_lib

shape_info = shape_info_lib.ShapesInfo(shape_curve_ids, bndry_verts, path_pairs)
print("shape_info.bndry_c_verts_start: ", shape_info.bndry_c_verts_start)
# print("9: ", shape_info.bndry_verts_flat[9])
dsplc_bndry_calc = dsplc_bndry_lib.DsplcCalc(shape_info)
dsplc_bndry_calc.get_pair_bndry_dsplcmnt()

print("dsplc_bndry_calc.pair_bndry_dsplcmnt: \n")
for i in range(len(dsplc_bndry_calc.pair_bndry_dsplcmnt)):
    print(path_pairs[i], ": ", dsplc_bndry_calc.pair_bndry_dsplcmnt[i])
dsplc_bndry_calc.get_var_ids()
print("dsplc_bndry_calc.variable_ids: ", dsplc_bndry_calc.variable_ids)
dsplc_bndry_calc.group_vars()

print("dsplc_bndry_calc.var_groups: ", dsplc_bndry_calc.var_groups)
dsplc_bndry_calc.buildGroupConnGraph()
# order_calc = order_calc_lib.OrderCalc(dsplc_bndry_calc.group_graph, dsplc_bndry_calc.group_closed)
# order_calc.determine_order()
dsplc_bndry_calc.getOrder()
# print("dsplc_bndry_calc.group_graph: ", dsplc_bndry_calc.group_graph)

print("group_bndry_map: ", dsplc_bndry_calc.order_calc.group_bndry_map)
print("ordered_group: ", dsplc_bndry_calc.order_calc.ordered_group)
print("dsplc_bndry_calc.group_closed: ", dsplc_bndry_calc.group_closed)

dsplc_bndry_calc.calcVariables()

print("dsplc_bndry_calc.pair_bndry_dsplcmnt: \n")
for dsplcmnt in dsplc_bndry_calc.pair_bndry_dsplcmnt:
    print(dsplcmnt)

dsplc_bndry_calc.itrpRestOfCurves()

write_lib.write_pts_out(shape_info, "bndry_pts.txt")
write_lib.write_uvs_out(shape_info, "bndry_pts_uv.txt")
write_lib.write_info_out(shape_info, dsplc_bndry_calc, "info.txt", seg_len, shape_verts_num)

# import matplotlib.pyplot as plt
# plt.figure(figsize=(20,10))
# plt.plot(shape_info.bndry_verts_flat[:,0], shape_info.bndry_verts_flat[:,1], 'o')
# plt.plot(shape_info.bndry_vert_uv[:,0], shape_info.bndry_vert_uv[:,1], 'o')
# # for i in range(len(shape_info.bndry_c_verts_start)-1):
# #     plt.plot(shape_info.bndry_vert_uv[shape_info.bndry_c_verts_start[i]:shape_info.bndry_c_verts_start[i+1],0], -shape_info.bndry_vert_uv[shape_info.bndry_c_verts_start[i]:shape_info.bndry_c_verts_start[i+1],1], 'o')
# #     plt.plot(shape_info.bndry_vert_uv[shape_info.bndry_c_verts_start[i], 0], -shape_info.bndry_vert_uv[shape_info.bndry_c_verts_start[i],1], 'o')
# # v=87
# # plt.plot(shape_info.bndry_vert_uv[v,0], -shape_info.bndry_vert_uv[v,1], 'o')
# plt.axis('scaled')
# plt.show()

# print("shape_info.path_pair: ", shape_info.path_pairs)
# print("group_bndry_map: ", order_calc.group_bndry_map)

# for i in range(len(dsplc_bndry_calc.var_groups)):
#     print(i, ": ", dsplc_bndry_calc.var_groups[i])
# print("dsplc_bndry_calc.var_groups: ", dsplc_bndry_calc.var_groups)
# print("dsplc_bndry_calc.group_closed: ", dsplc_bndry_calc.group_closed)
# for r in range(len(dsplc_bndry_calc.group_graph)):
# print(dsplc_bndry_calc.group_graph[0])
