import gurobipy as gp
from gurobipy import GRB

angles = []
angle_diffs_constr_dic = {}
shape_curve_ids = []

f = open("polygno_out.txt", "r")
lines = f.readlines()

cnt = 1
while not lines[cnt].startswith("angle_diffs_constr_dic"):
    segs = lines[cnt].split()
    angles.append(float(segs[0]))
    cnt += 1
cnt += 1
while not lines[cnt].startswith("shape_curve_ids"):
    segs = lines[cnt].split()
    angle_diffs_constr_dic[(int(segs[0]), int(segs[1]))] = float(segs[2])
    cnt += 1
cnt += 1
while cnt < len(lines):
    segs = lines[cnt].split()
    cids = []
    for seg in segs:
        cids.append(int(seg))
    shape_curve_ids.append(cids)
    cnt += 1

# print("angle_diffs_dic: ", angle_diffs_dic)
# print("angle_diffs_constr_dic: ", angle_diffs_constr_dic)

class Opt:
    def __init__(self, angle_diffs_constr_dic, angles, shape_curve_ids):
        self.angles = angles
        self.shape_curve_ids = shape_curve_ids
        self.angle_diffs_constr_dic = angle_diffs_constr_dic
        self.m = gp.Model("model")

    def createVars(self):
        self.a = [None] * len(self.shape_curve_ids)
        self.diff = [None] * len(self.angle_diffs_constr_dic)
        self.d = [None] * len(self.angle_diffs_constr_dic)

        for i in range(len(self.a)):
            self.a[i] = self.m.addVar(name="a_"+str(i), lb = 0, ub = 360)
        for i in range(len(self.diff)):
            self.diff[i] = self.m.addVar(name="diff_"+str(i), lb = -360, ub = 360)
        for i in range(len(self.d)):
            self.d[i] = self.m.addVar(name="k_"+str(i), vtype = GRB.INTEGER, lb = -4, ub = 4)

    def getShapeId(self, cid):
        for shape_id in range(len(self.shape_curve_ids)):
            if cid in self.shape_curve_ids[shape_id]:
                return shape_id

    def addConstraints(self):
        cnt = 0
        d_cnt = 0
        for key in self.angle_diffs_constr_dic:
            shape_id_0 = self.getShapeId(key[0])
            shape_id_1 = self.getShapeId(key[1])
            self.m.addConstr( self.diff[cnt] == self.a[shape_id_0] + self.angles[key[0]] - self.a[shape_id_1] - self.angles[key[1]] - self.d[d_cnt] * 90)
            cnt += 1
            d_cnt += 1

    def solve(self):
        objective = gp.QuadExpr()
        for i in range(len(self.diff)):
            objective += self.diff[i] * self.diff[i]

        # Set objective: maximize x
        self.m.setObjective(objective, GRB.MINIMIZE)

        try:
            # self.m.write("m.lp")
            # m.feasRelaxS(0, True, False, True)
            # m.write("feasopt1.lp")
            self.m.optimize()
        except gp.GurobiError:
            print("Optimize failed due to non-convexity")

        f = open("optmized_angles.txt", "w")
        for i in range(len(self.a)):
            print("a_", i, self.a[i].x)
            f.write(str(self.a[i].x) + "\n")
        f.close()
        # for v in self.m.getVars():
        #     print('%s %g' % (v.varName, v.x))
        # f = open("optmized_rels.txt", "w")
        # for i in range(len(self.d)):
        #     print("d_", i, self.d[i].x)
        #     f.write(str(self.d[i].x) + "\n")
        # f.close()

opt = Opt(angle_diffs_constr_dic, angles, shape_curve_ids)
opt.createVars()
opt.addConstraints()
opt.solve()
