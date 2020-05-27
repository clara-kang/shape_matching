import numpy as np

class ShapesInfo:

    def flatten_bndry_verts(self):
        self.bndry_verts_flat = []
        self.bndry_c_verts_start = np.zeros(len(self.bndry_verts))
        cnt = 0
        for curve_id in range(len(self.bndry_verts)):
            curve_verts = self.bndry_verts[curve_id]
            self.bndry_verts_flat.extend(curve_verts)
            self.bndry_c_verts_start[curve_id] = cnt
            cnt += len(curve_verts)

        self.bndry_verts_flat = np.array(self.bndry_verts_flat)
        self.bndry_c_verts_start = self.bndry_c_verts_start.astype(int)

    def  __init__(self, shape_curve_ids, bndry_verts, path_pairs):
        self.shape_curve_ids = shape_curve_ids
        self.bndry_verts = bndry_verts
        self.path_pairs = path_pairs
        self.flatten_bndry_verts()
        self.bndry_vert_uv = bndry_vert_uv = np.copy(self.bndry_verts_flat)

    def getNextCurveInShape(self, c_id):
        for shape_id in range(len(self.shape_curve_ids)):
            # find the shape that contains c_id
            if c_id in self.shape_curve_ids[shape_id]:
                c_in_shape_idx = self.shape_curve_ids[shape_id].index(c_id)
                return self.shape_curve_ids[shape_id][(c_in_shape_idx + 1) % len(self.shape_curve_ids[shape_id])]

    def getPrevCurveInShape(self, c_id):
        for shape_id in range(len(self.shape_curve_ids)):
            # find the shape that contains c_id
            if c_id in self.shape_curve_ids[shape_id]:
                c_in_shape_idx = self.shape_curve_ids[shape_id].index(c_id)
                return self.shape_curve_ids[shape_id][(c_in_shape_idx - 1 + len(self.shape_curve_ids[shape_id])) % len(self.shape_curve_ids[shape_id])]

    def getPairCurveVids(self, pair_info):
        c_id1 = pair_info[0]
        c_id2 = pair_info[1]
        curve_vid_1 = np.arange(self.bndry_c_verts_start[c_id1], self.bndry_c_verts_start[c_id1] + len(self.bndry_verts[c_id1]))
        curve_vid_2 = np.arange(self.bndry_c_verts_start[c_id2], self.bndry_c_verts_start[c_id2] + len(self.bndry_verts[c_id2]))
        curve_vid_1 = np.append(curve_vid_1, self.bndry_c_verts_start[self.getNextCurveInShape(c_id1)]).astype(int)
        curve_vid_2 = np.append(curve_vid_2, self.bndry_c_verts_start[self.getNextCurveInShape(c_id2)]).astype(int)

        if not (pair_info[2]):
            curve_vid_2 = np.flip(curve_vid_2)

        return (curve_vid_1, curve_vid_2)

    def getCurveEnd(self, cid):
        next_cid = self.getNextCurveInShape(cid)
        return self.bndry_c_verts_start[next_cid]

    def getEPs(self, cid):
        return [self.bndry_c_verts_start[cid], self.getCurveEnd(cid)]

    def get_ep_pair(self, pair_id):
        def swap(arr):
            return [arr[1], arr[0]]

        path_pair = self.path_pairs[pair_id]

        eps = self.getEPs(path_pair[0])
        ep_pairs = self.getEPs(path_pair[1])
        if not (path_pair[2]):
            ep_pairs = swap(ep_pairs)

        return (eps, ep_pairs)
