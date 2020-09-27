import numpy as np
import math

class HE:
    def __init__(self, vert, face, he_next=None, pair=None):
        self.vert = vert
        self.pair = pair
        self.face = face
        self.he_next = he_next



def map(faces, vertices, bndry_vals):

    def getPrevHeId(he_id):
        nonlocal HEs
        return HEs[HEs[he_id].he_next].he_next


    def getNbs(v_id) :
        nonlocal HEs, HE_verts
        neighbors = []
        he = HE_verts[v_id]
        he = he.he_next
        first_he = he

        first_nb = he.vert
        neighbors.append(first_nb)

        if he.pair != None:
            he = he.pair.he_next
            on_bndry = False

            while he.vert != first_nb:
                neighbors.append(he.vert)
                if he.pair == None:
                    on_bndry = True
                    break
                he = he.pair.he_next
        else:
            on_bndry = True

        if on_bndry:
            other_side_nbs = []
            he =  first_he.he_next
            other_side_nbs.append(he.vert)
            he = he.he_next
            if he.pair != None:
                he = he.pair.he_next
                cnt = 0
                while he.vert != other_side_nbs[0]:
                    other_side_nbs.append(he.vert)
                    he = he.he_next
                    if he.pair == None:
                        cnt = 100
                        break
                    he = he.pair.he_next
                    cnt += 1
            other_side_nbs.reverse()
            other_side_nbs.extend(neighbors)
            return other_side_nbs
        return neighbors


    # get weight for half edge
    def getHEWeight(he):
        nonlocal vertices, HEs
        vert_3 = np.array(vertices[he.he_next.vert])
        vert_1 = np.array(vertices[he.vert])
        vert_2 = np.array(vertices[he.he_next.he_next.vert])
        v1 = vert_1 - vert_3
        v2 = vert_2 - vert_3
    #     print("v1, ", v1, ", v2, ", v2)
        cos_angle = np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2))
        sign = np.cross(v1, v2)
    #     print("sign: ", sign)
        # coordinate system has z axis going into the screen
        if sign > 0:
            cos_angle = -cos_angle
        sin_angle = math.sqrt(1 - cos_angle * cos_angle)
        return cos_angle / sin_angle


    def getEdgeWeight(start_v, end_v):
        nonlocal pairs
        weight = 0
        if (start_v, end_v) in pairs:
            weight += getHEWeight(pairs[(start_v, end_v)])
        if (end_v, start_v) in pairs:
            weight += getHEWeight(pairs[(end_v, start_v)])
        return weight


    HEs = [None] * (len(faces) * 3)
    HE_verts = [None] * len(vertices)
    pairs = {}
    bndry_n = len(bndry_vals)
    N = len(vertices)

    # build HE structure
    for face_id in range (0, len(faces)):
        face = faces[face_id]
        face_HEs = []
        for i in range (0, 3):
            v_id = face[i]
            v_next_id = face[(i+1)%3]
            he = HE(v_next_id, face_id)
            face_HEs.append(he)
            HEs[int(3 * face_id + i)] = he
            pairs[(v_id, v_next_id)] = he
            HE_verts[v_next_id] = he
        for i in range(0, 3):
            face_HEs[i].he_next = face_HEs[(i+1)%3]

    # find pairs
    for he_id in range(len(HEs)):
        # get start, end vertex for he
        end_v = HEs[he_id].vert
        start_v = HEs[he_id].he_next.he_next.vert
        if (end_v, start_v) in pairs:
            HEs[he_id].pair = pairs[(end_v, start_v)]

    M = np.zeros([N - bndry_n, N - bndry_n])
    rhs = np.zeros(N - bndry_n)

    # build M
    for v_idx in range (bndry_n, N):
        mtrx_idx = v_idx - bndry_n
        nbs = getNbs(v_idx)
        nb_weights = np.zeros(len(nbs))
        for nb_idx in range(len(nbs)):
            nb = nbs[nb_idx]
            nb_weights[nb_idx] = getEdgeWeight(v_idx, nb)
        nb_weights = nb_weights / np.sum(nb_weights)
        M[mtrx_idx][mtrx_idx] = 1

        for nb_idx in range(len(nbs)):
            nb = nbs[nb_idx]
            # on boundary
            if nb < bndry_n:
                rhs[mtrx_idx] += bndry_vals[nb] * nb_weights[nb_idx]
            else:
                M[mtrx_idx][nb - bndry_n] -= nb_weights[nb_idx]

    sol = np.linalg.solve(M, rhs)
    vals = np.concatenate((bndry_vals, sol))
    return vals
