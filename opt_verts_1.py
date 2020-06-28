vert_locs = []
shape_verts = []
pair_rels = []
f = open("poly_verts_rot.txt", "r")
for line in f.readlines():
    segs = line.split()
    vert_locs.append([float(segs[0]), float(segs[1])])
f.close()

f =  open("polygon_out_2.txt", "r")
lines = f.readlines()
cnt = 1
while not lines[cnt].startswith("pair_rels:"):
    segs = lines[cnt].split()
    vids = []
    for seg in segs:
        vids.append(int(seg))
    shape_verts.append(vids)
    cnt += 1
cnt += 1
while cnt < len(lines):
    segs = lines[cnt].split()
    pair_rels.append([int(segs[0]), int(segs[1]), int(segs[2]), int(segs[3]), int(segs[4])])
    cnt += 1

import gurobipy as gp
from gurobipy import GRB

class OptVerts:
    def __init__(self, vert_locs, shape_verts, pair_rels):
        self.vert_locs = vert_locs
        self.pair_rels = pair_rels
        self.shape_verts = shape_verts
        # self.px = [None] * len(vert_locs)
        # self.py = [None] * len(vert_locs)
        # self.diff_x = [None] * len(vert_locs)
        # self.diff_y = [None] * len(vert_locs)
        # self.std_diff_x = [None] * len(vert_locs)
        # self.std_diff_y = [None] * len(vert_locs)
        # self.abs_std_diff_x = [None] * len(vert_locs)
        # self.abs_std_diff_y = [None] * len(vert_locs)
        # self.ds_x = [None] * len(pair_rels)
        # self.ds_y = [None] * len(pair_rels)

    def createVars(self):
        self.m = gp.Model("model")
        self.px = self.m.addVars(len(self.vert_locs), lb = -GRB.INFINITY, name="px")
        self.py = self.m.addVars(len(self.vert_locs), lb = -GRB.INFINITY, name="py")
        self.diff_x = self.m.addVars(len(self.vert_locs), lb = -GRB.INFINITY, name="diff_x")
        self.diff_y = self.m.addVars(len(self.vert_locs), lb = -GRB.INFINITY, name="diff_y")
        self.std_diff_x = self.m.addVars(len(self.vert_locs), lb = -GRB.INFINITY, name="std_diff_x")
        self.std_diff_y = self.m.addVars(len(self.vert_locs), lb = -GRB.INFINITY, name="std_diff_y")
        self.abs_std_diff_x = self.m.addVars(len(self.vert_locs), name="abs_std_diff_x")
        self.abs_std_diff_y = self.m.addVars(len(self.vert_locs), name="abs_std_diff_y")

        self.ds_x = self.m.addVars(len(pair_rels)-1, lb = -GRB.INFINITY, vtype = GRB.INTEGER, name="ds_x")
        self.ds_y = self.m.addVars(len(pair_rels)-1, lb = -GRB.INFINITY, vtype = GRB.INTEGER, name="ds_y")

    def createConstraints(self):
        for i in range(len(self.vert_locs)):
            self.m.addConstr(self.diff_x[i] == self.px[i] - self.vert_locs[i][0])
            self.m.addConstr(self.diff_y[i] == self.py[i] - self.vert_locs[i][1])
        for i in range(len(self.shape_verts)):
            for vid in self.shape_verts[i]:
                self.m.addConstr(self.diff_x[vid] - sum(self.diff_x[vid_1] for vid_1 in self.shape_verts[i]) / len(self.shape_verts[i]) == self.std_diff_x[vid])
                self.m.addConstr(self.abs_std_diff_x[vid] == gp.abs_(self.std_diff_x[vid]))
                self.m.addConstr(self.diff_y[vid] - sum(self.diff_y[vid_1] for vid_1 in self.shape_verts[i]) / len(self.shape_verts[i]) == self.std_diff_y[vid])
                self.m.addConstr(self.abs_std_diff_y[vid] == gp.abs_(self.std_diff_y[vid]))
        for i in range(len(pair_rels)-1):
            p_rel = pair_rels[i]
            segs = [[p_rel[0], p_rel[1]], [p_rel[2], p_rel[3]]]
            for se in range(2):
                # (x, y) -> (x, y)
                if p_rel[4] == 0:
                    self.m.addConstr(self.px[segs[0][se]] - self.px[segs[1][se]] == self.ds_x[i])
                    self.m.addConstr(self.py[segs[0][se]] - self.py[segs[1][se]] == self.ds_y[i])
                # (x, y) -> (-y, x)
                elif p_rel[4] == 1:
                    self.m.addConstr(-self.py[segs[0][se]] - self.px[segs[1][se]] == self.ds_x[i])
                    self.m.addConstr(self.px[segs[0][se]] - self.py[segs[1][se]] == self.ds_y[i])
                # (x, y) -> (y, -x)
                elif p_rel[4] == 2:
                    self.m.addConstr(self.py[segs[0][se]] - self.px[segs[1][se]] == self.ds_x[i])
                    self.m.addConstr(-self.px[segs[0][se]] - self.py[segs[1][se]] == self.ds_y[i])
                # (x, y) -> (x, y)
                elif p_rel[4] == 3:
                    self.m.addConstr(-self.px[segs[0][se]] - self.px[segs[1][se]] == self.ds_x[i])
                    self.m.addConstr(-self.py[segs[0][se]] - self.py[segs[1][se]] == self.ds_y[i])

    def addObjective(self):
        # self.objective = gp.LinExpr()
        # for i in range(len(self.vert_locs)):
        #     self.objective += self.abs_std_diff_x[i]
        #     self.objective += self.abs_std_diff_y[i]
        #
        # self.m.setObjective(self.objective, GRB.MINIMIZE)
        self.m.setObjective(sum(self.abs_std_diff_x[i] for i in range(len(self.abs_std_diff_x))) + \
            sum(self.abs_std_diff_y[i] for i in range(len(self.abs_std_diff_y))), GRB.MINIMIZE)

    def solve(self):
        try:
            # m.params.NonConvex = 2
            self.m.write("opt_v.lp")
            # self.m.feasRelaxS(0, True, False, True)
            # self.m.write("opt_verts_relaxed.lp")
            self.m.optimize()
        except gp.GurobiError:
            print("Optimize failed")

        #
        f = open("opt_pts.txt", "w")
        for i in range(len(self.vert_locs)):
            print("x_", i, ": ", self.px[i].x)
            print("y_", i, ": ", self.py[i].x)
            f.write(str(self.px[i].x) + " " + str(self.py[i].x) + "\n")
        f.close()

        for v in self.m.getVars():
            print('%s %g' % (v.varName, v.x))

opt = OptVerts(vert_locs, shape_verts, pair_rels)
opt.createVars()
opt.createConstraints()
opt.addObjective()
opt.solve()
