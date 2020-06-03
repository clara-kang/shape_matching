
def write_pts_out(shape_info, file_name):
    pts_file = open(file_name, "w")

    for curve_id in range(len(shape_info.bndry_verts)):
        for v in shape_info.bndry_verts[curve_id]:
            pts_file.write(str(v[0]) + " " + str(v[1]) + "\n")

    pts_file.close()

def write_uvs_out(shape_info, file_name):
    pts_uv_file = open(file_name, "w")

    for uv in shape_info.bndry_vert_uv:
        pts_uv_file.write(str(uv[0]) + " " + str(uv[1]) + "\n")

    pts_uv_file.close()

def write_info_out(shape_info, dsplc_bndry_calc, file_name, seg_len, shape_verts_num):
    info_file = open(file_name, "w")

    info_file.write("seg_len="+str(seg_len)+"\n")
    info_file.write("bndry_pts_cnt=")

    shape_verts_num_str = ""
    for num in shape_verts_num:
        shape_verts_num_str += str(num) + ", "
    shape_verts_num_str = shape_verts_num_str[:-2]
    info_file.write(shape_verts_num_str)

    info_file.write("\nbndry_c_verts_start=")
    c_v_start_str = ""
    for c_v_start in shape_info.bndry_c_verts_start:
        c_v_start_str += str(c_v_start) + ", "
    c_v_start_str = c_v_start_str[:-2]
    info_file.write(c_v_start_str)

    info_file.write("\npath_pair=")
    pair_str = ""
    for pair in shape_info.path_pairs:
        pair_str += "(" + str(pair[0]) + ", " + str(pair[1]) + ", " + str(int(pair[2])) + ")" + ", "
    pair_str = pair_str[:-2]
    info_file.write(pair_str)

    info_file.write("\npair_bndry_dsplcmnt=")
    for p_dsplcmnt in dsplc_bndry_calc.pair_bndry_dsplcmnt:
        func_type = p_dsplcmnt[1]
        if func_type == 2:
            func_type = 1
        elif func_type == 1:
            func_type = 2
        info_file.write("[" + str(p_dsplcmnt[0][0]) + ", " + str(p_dsplcmnt[0][1]) + ", " + str(func_type) + "]" + ", ")

    info_file.write("\nshape_curve_ids=")
    for sc_id in shape_info.shape_curve_ids:
        info_file.write(str(sc_id) + ", ")

    info_file.close()
