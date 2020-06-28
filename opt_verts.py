vert_locs = []
shape_verts = []
pair_rels = []
f = open("poly_verts.txt", "r")
for line in f.readlines():
    segs = line.split()
    vert_locs.append([float(segs[0]), float(segs[1])])
f.close()

f =  open("polygon_out_2.txt", "r")
lines = f.readlines()
cnt = 2
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
        self.px = [None] * len(vert_locs)
        self.py = [None] * len(vert_locs)
        self.diff_x = [None] * len(vert_locs)
        self.diff_y = [None] * len(vert_locs)
        self.std_diff_x = [None] * len(vert_locs)
        self.std_diff_y = [None] * len(vert_locs)
        self.abs_std_diff_x = [None] * len(vert_locs)
        self.abs_std_diff_y = [None] * len(vert_locs)
        self.ds_x = [None] * len(pair_rels)
        self.ds_y = [None] * len(pair_rels)

    def createVars(self):
        self.m = gp.Model("model")
        for i in range(len(self.vert_locs)):
            self.px[i] = self.m.addVar(name="px_"+str(i))
            self.py[i] = self.m.addVar(name="py_"+str(i))
            self.diff_x[i] = self.m.addVar(name="diff_x_"+str(i), lb = -GRB.INFINITY)
            self.diff_y[i] = self.m.addVar(name="diff_y_"+str(i), lb = -GRB.INFINITY)
            self.std_diff_x[i] = self.m.addVar(name="std_diff_x_"+str(i), lb = -GRB.INFINITY)
            self.std_diff_y[i] = self.m.addVar(name="std_diff_y_"+str(i), lb = -GRB.INFINITY)
            self.abs_std_diff_x[i] = self.m.addVar(name="abs_std_diff_x_"+str(i))
            self.abs_std_diff_y[i] = self.m.addVar(name="abs_std_diff_y_"+str(i))
        for i in range(len(pair_rels)):
            self.ds_x[i] = self.m.addVar(name="ds_x"+str(i), vtype = GRB.INTEGER, lb = -GRB.INFINITY)
            self.ds_y[i] = self.m.addVar(name="ds_y"+str(i), vtype = GRB.INTEGER, lb = -GRB.INFINITY)

    def createConstraints(self):
        for i in range(len(self.vert_locs)):
            self.m.addConstr(self.diff_x[i] == self.px[i] - self.vert_locs[i][0])
            self.m.addConstr(self.diff_y[i] == self.py[i] - self.vert_locs[i][1])
        for i in range(len(self.shape_verts)):
            for vid in self.shape_verts[i]:
                # self.m.addConstr(self.diff_x[vid] - self.diff_x.sum(vid_1, '*') / len(self.shape_verts[i]) == self.std_diff_x[vid] for vid_1 in self.shape_verts[i])
                self.m.addConstr(self.diff_x.sum(vid_1, '*') == 0 for vid_1 in self.shape_verts[i])
                self.m.addConstr(self.abs_std_diff_x[vid] == gp.abs_(self.std_diff_x[vid]))
                self.m.addConstr(self.diff_y[vid] - self.diff_y.sum(vid_1, '*') / len(self.shape_verts[i]) == self.std_diff_y[vid] for vid_1 in self.shape_verts[i])
                self.m.addConstr(self.abs_std_diff_y[vid] == gp.abs_(self.std_diff_y[vid]))
        for i in range(len(pair_rels)):
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
            self.objective = gp.LinExpr()
            for i in range(self.vert_locs):
                self.objective += self.abs_std_diff_x[i]
                self.objective += self.abs_std_diff_y[i]

            self.m.setObjective(self.objective, GRB.MINIMIZE)

        def solve(self):
            try:
                # m.params.NonConvex = 2
                m.write("opt_verts.lp")
                m.feasRelaxS(0, True, False, True)
                m.write("opt_verts_relaxed.lp")
                m.optimize()
            except gp.GurobiError:
                print("Optimize failed")

            #
            f = open("opt_pts.txt", "w")
            for i in range(8):
                print("x_", i, ": ", px[i].x)
                print("y_", i, ": ", py[i].x)
                f.write(str(px[i].x) + " " + str(py[i].x) + "\n")
            f.close()

            for v in m.getVars():
                print('%s %g' % (v.varName, v.x))

opt = OptVerts(vert_locs, shape_verts, pair_rels)
opt.createVars()
opt.createConstraints()
opt.solve()
