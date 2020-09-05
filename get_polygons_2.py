
import numpy as np
import math
from random import randrange

# V2 - V1
def getSignedAngle(v1, v2):
    cos_angle = np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = min(cos_angle, 1)
    cos_angle = max(cos_angle, -1)
    angle = math.acos(cos_angle)
    sign = np.cross(v1, v2)
    if sign < 0:
        angle = -angle
        # angle += np.pi * 2
    return angle

def radToDegree(rad):
    return 180 * rad / np.pi

def degreeToRad(degree):
    return degree * np.pi / 180

def flatten_arr(arr):
    flat_arr = []
    for i in range(len(arr)):
        flat_arr.extend(arr[i])
    return flat_arr

def angleNormalize(angle):
    angle = angle % 360
    if angle < 0:
        angle += 360
    return angle

def rayIntrsect(p1, r1, p2, r2):
    M = np.array([[r1[0], -r2[0]], [r1[1], -r2[1]]])
    B = p2 - p1
    ts = np.linalg.solve(M, B) # [t1, t2]
    return ts[0]

def rotate(rot_angle, v):
    rot_angle_rad = degreeToRad(rot_angle)
    rot_matrix = np.array([[math.cos(rot_angle_rad), -math.sin(rot_angle_rad)], [math.sin(rot_angle_rad), math.cos(rot_angle_rad)]])
    rot_v = np.matmul(rot_matrix, v)
    return rot_v

class Polygons:
    def __init__(self, bndry_verts, shape_curve_ids, path_pairs):
        self.shape_curve_ids = shape_curve_ids
        self.bndry_verts = bndry_verts
        self.bndry_verts_uv = [None] * len(self.bndry_verts)
        self.path_pairs = path_pairs
        # print("path_pairs: ", path_pairs)
        self.pair_func_id = [None] * len(path_pairs)
        self.poly_verts = []
        self.poly_vecs = []
        self.poly_norms = []
        self.shape_ccw = [False] * len(self.shape_curve_ids)
        self.poly_centers = [None] * len(self.shape_curve_ids)
        self.poly_vecs_flat = []
        self.angles = []
        self.angles_flat = []
        self.angle_diffs_constr_dic = {}

    def getCid(self, shape_id, which_id):
        cnt = 0
        for i in range(shape_id):
            cnt += len(self.shape_curve_ids[i])
        cnt += which_id
        return cnt

    def getWhichId(cid):
        cnt = 0
        shape_id = 0
        while cnt + self.shape_curve_ids[shape_id] < cid:
            cnt = cnt + self.shape_curve_ids[shape_id]
            shape_id += 1
        return shape_id, cid - cnt

    def getPolyVerts(self):
        for i in range(len(self.shape_curve_ids)):
            self.poly_verts.append([])
            for cid in self.shape_curve_ids[i]:
                self.poly_verts[-1].append(np.array(self.bndry_verts[cid][0]))

    def applyMultiplier(self, multiplier):
        for i in range(len(self.bndry_verts)):
            for j in range(len(self.bndry_verts[i])):
                self.bndry_verts[i][j] *= multiplier

    def writeVertsOut(self, file_name):
        f = open(file_name, "w")
        for cid in range(len(self.bndry_verts)):
            v = np.array(self.bndry_verts[cid][0])
            f.write(str(v[0]) + " " + str(v[1]) + "\n")
        f.close()

    def getPolyVecs(self):
        self.poly_vecs = []
        self.poly_vecs_flat = [None] * len(self.bndry_verts)
        for i in range(len(self.poly_verts)):
            self.poly_vecs.append([])
            for j in range(len(self.poly_verts[i])):
                vec = self.poly_verts[i][(j+1) % len(self.poly_verts[i])] - self.poly_verts[i][j]
                self.poly_vecs[-1].append(vec)
        for cid in range(len(self.bndry_verts)):
            seg = self.getStartEndV(cid)
            self.poly_vecs_flat[cid] = self.bndry_verts[seg[1]][0] - self.bndry_verts[seg[0]][0]

    def getCenters(self):
        for i in range(len(self.poly_vecs)):
            self.poly_centers[i] = np.zeros(2)
            for j in range(len(self.poly_verts[i])):
                self.poly_centers[i] += self.poly_verts[i][j]
            self.poly_centers[i] /= len(self.poly_verts[i])

        print("len(self.poly_centers): ", len(self.poly_centers))
        for i in range(len(self.poly_centers)):
            print(self.poly_centers[i])
        print("self.poly_centers: ", self.poly_centers)

    def getAngles(self):
        self.angles = []
        for i in range(len(self.poly_vecs)):
            self.angles.append([])
            for j in range(len(self.poly_vecs[i])):
                self.angles[-1].append(angleNormalize(radToDegree(getSignedAngle(np.array([1.0, 0.0]), self.poly_vecs[i][j]))))
        for i in range(len(self.poly_vecs)):
            for j in range(len(self.poly_vecs[i])):
                cid1 = self.getCid(i, j)
                cid2 = self.getCid(i, (j+1)%len(self.poly_vecs[i]))
        self.angles_flat = flatten_arr(self.angles)

    def getPairAngleDiff(self):
        self.angle_diffs_constr_dic = {}
        for i in range(len(self.path_pairs)):
            ppair = self.path_pairs[i]
            vec1 = self.poly_vecs_flat[ppair[0]]
            vec2 = self.poly_vecs_flat[ppair[1]]
            if ppair[2] == False:
                vec2 = - vec2
            angle_diff = radToDegree(getSignedAngle(vec1, vec2))
            positive = angle_diff > 0
            rounded_angle = round(abs(angle_diff) / 90) * 90
            if not positive:
                rounded_angle = -rounded_angle
            self.angle_diffs_constr_dic[(ppair[1], ppair[0])] = rounded_angle

            if rounded_angle == 0:
                self.pair_func_id[i] = 0
            elif rounded_angle == 90:
                self.pair_func_id[i] = 1
            elif rounded_angle == -90:
                self.pair_func_id[i] = 2
            elif rounded_angle == 180 or rounded_angle == -180:
                self.pair_func_id[i] = 3
            else:
                print("rounded_angle: ", rounded_angle)

    def writeOut(self):
        f = open("polygno_out.txt", "w")
        f.write("angles\n")
        for angle in self.angles_flat:
            f.write(str(angle) + "\n")
        f.write("angle_diffs_constr_dic\n")
        for key in self.angle_diffs_constr_dic:
            f.write(str(key[0]) + " " + str(key[1]) + " " + str(self.angle_diffs_constr_dic[key]) + "\n")
        f.write("shape_curve_ids\n")
        for cids in self.shape_curve_ids:
            for cid in cids:
                f.write(str(cid) + " ")
            f.write("\n")

    def getNextCurveInShape(self, c_id):
        for shape_id in range(len(self.shape_curve_ids)):
            # find the shape that contains c_id
            if c_id in self.shape_curve_ids[shape_id]:
                c_in_shape_idx = self.shape_curve_ids[shape_id].index(c_id)
                return self.shape_curve_ids[shape_id][(c_in_shape_idx + 1) % len(self.shape_curve_ids[shape_id])]

    def getStartEndV(self, cid):
        next_cid = self.getNextCurveInShape(cid)
        return (cid, next_cid)

    def writeOut2(self):
        f = open("polygon_out_2.txt", "w")
        f.write("shape_curve_ids: \n")
        for i in range(len(self.shape_curve_ids)):
            for vid in self.shape_curve_ids[i]:
                f.write(str(vid) + " ")
            f.write("\n")
        f.write("pair_rels: \n")
        for i in range(len(self.path_pairs)):
            pair = self.path_pairs[i]
            seg0 = self.getStartEndV(pair[0])
            seg1 = self.getStartEndV(pair[1])
            if pair[2] == False:
                seg1 = (seg1[1], seg1[0])
            f.write(str(seg0[0]) + " " + str(seg0[1]) + " " + str(seg1[0]) + " " + str(seg1[1]) + " " + str(self.pair_func_id[i]) + "\n")


    def readInRotation(self):
        f = open("optmized_angles.txt", "r")
        self.shape_rotation = []
        cnt = 0
        for line in f.readlines():
            angle = float(line)
            angle = angleNormalize(angle)
            self.shape_rotation.append(angle)
        f.close()
        print("self.shape_rotation: ", self.shape_rotation)

    def readInVerts(self):
        f = open("opt_pts.txt", "r")
        self.optimized_verts = []
        for line in f.readlines():
            segs = line.split()
            self.optimized_verts.append(np.array([float(segs[0]), float(segs[1])]))
        f.close()

    def cidInPair(self, cid):
        for ppair in self.path_pairs:
            if cid == ppair[0] or cid == ppair[1]:
                return True
        return False

    def getCurveLen(self, cid):
        next_cid = self.getNextCurveInShape(cid)
        end_v = self.bndry_verts[next_cid][0]
        length = 0
        for i in range(1, len(self.bndry_verts[cid])):
            length += np.linalg.norm(np.array(self.bndry_verts[cid][i]) - np.array(self.bndry_verts[cid][i-1]))
        length += np.linalg.norm(end_v - np.array(self.bndry_verts[cid][-1]))
        return length

    def getCurvePartLen(self, cid, stop_vid):
        length = 0
        for i in range(1, stop_vid+1):
            length += np.linalg.norm(np.array(self.bndry_verts[cid][i]) - np.array(self.bndry_verts[cid][i-1]))
        return length

    def intrpCurve(self):
        for cid in range(len(self.bndry_verts)):
            v_num = len(self.bndry_verts[cid])
            curve_verts = [None] * v_num
            next_cid = self.getNextCurveInShape(cid)
            start_v = self.optimized_verts[cid]
            end_v = self.optimized_verts[next_cid]
            curve_verts[0] = start_v
            # check if cid in pair
            curve_len = self.getCurveLen(cid)
            if self.cidInPair(cid):
                # print("in pair")
                for i in range(1, v_num):
                    # t = self.getCurvePartLen(cid, i) / curve_len
                    t = i / v_num
                    v = (1 - t) * start_v + t * end_v
                    curve_verts[i] = v
            # curve is free
            else:
                # print("not in pair")
                start_displcmnt = start_v - self.bndry_verts[cid][0]
                end_displcmnt = end_v - self.bndry_verts[next_cid][0]
                for i in range(1, v_num):
                    t = i / v_num
                    v = (1 - t) * start_displcmnt + t * end_displcmnt + self.bndry_verts[cid][i]
                    curve_verts[i] = v
            self.bndry_verts_uv[cid] = curve_verts

    def rotateShapes(self):
        print("self.angles: ", self.angles)
        for shape_id in range(len(self.angles)):
            rot_angle = self.shape_rotation[shape_id]
            print("rot_angle: ", rot_angle)
            rot_angle_rad = degreeToRad(rot_angle)
            rot_matrix = np.array([[math.cos(rot_angle_rad), -math.sin(rot_angle_rad)], [math.sin(rot_angle_rad), math.cos(rot_angle_rad)]])

            for vid in range(len(self.angles[shape_id])):
                v = self.poly_verts[shape_id][vid]
                rot_v = np.matmul(rot_matrix, v)
                self.poly_verts[shape_id][vid] = rot_v

            # get new shape_center
            new_shape_center = np.zeros(2)
            for vid in range(len(self.angles[shape_id])):
                new_shape_center += self.poly_verts[shape_id][vid]
            new_shape_center /= len(self.angles[shape_id])

            # displace to old center
            disp = self.poly_centers[shape_id] - new_shape_center
            print("disp: ", disp)
            for vid in range(len(self.angles[shape_id])):
                self.poly_verts[shape_id][vid] += disp

            for cid in self.shape_curve_ids[shape_id]:
                for i in range(len(self.bndry_verts[cid])):
                    self.bndry_verts[cid][i] = np.matmul(rot_matrix, self.bndry_verts[cid][i])
                    self.bndry_verts[cid][i] += disp

    def getRandomShapeRotations(self):
        self.shape_rotation = [None] * len(self.angles)
        for i in range(len(self.shape_rotation)):
            self.shape_rotation[i] = 90 * randrange(4)

    def resetShapeRotations(self):
        self.shape_rotation = [0] * len(self.angles)

    def get_bndry_c_verts_start(self):
        cnt = 0
        bndry_c_verts_start = []
        for curve in self.bndry_verts:
            bndry_c_verts_start.append(cnt)
            cnt += len(curve)
        return bndry_c_verts_start
