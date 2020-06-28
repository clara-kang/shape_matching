import gurobipy as gp
from gurobipy import GRB

angle_num = 0
angle_diffs_dic = {}
angle_diffs_constr_dic = {}

f = open("polygno_out.txt", "r")
lines = f.readlines()
angle_num = int(lines[0])

cnt = 2
while not lines[cnt].startswith("angle_diffs_constr_dic"):
    segs = lines[cnt].split()
    angle_diffs_dic[(int(segs[0]), int(segs[1]))] = float(segs[2])
    cnt += 1
cnt += 1
while cnt < len(lines):
    segs = lines[cnt].split()
    angle_diffs_constr_dic[(int(segs[0]), int(segs[1]))] = float(segs[2])
    cnt += 1

# print("angle_diffs_dic: ", angle_diffs_dic)
# print("angle_diffs_constr_dic: ", angle_diffs_constr_dic)

class Opt:
    def __init__(self, var_num, angle_diffs_constr_dic, angle_diffs_dic):
        self.angle_diffs_dic = angle_diffs_dic
        self.var_num = var_num
        self.angle_diffs_constr_dic = angle_diffs_constr_dic
        self.m = gp.Model("model")

    def createVars(self):
        self.a = [None] * self.var_num
        self.k = [None] * ( len(self.angle_diffs_dic) + len(self.angle_diffs_constr_dic))
        self.diff = [None] * len(self.angle_diffs_dic)
        self.abs_diff = [None] * len(self.angle_diffs_dic)

        for i in range(self.var_num):
            self.a[i] = self.m.addVar(name="a_"+str(i), lb = -360, ub = 360)
        for i in range(len(self.angle_diffs_dic)):
            self.diff[i] = self.m.addVar(name="diff_"+str(i), lb = -360, ub = 360)
            self.abs_diff[i] = self.m.addVar(name="abs_diff_"+str(i), ub = 360)
        for i in range(len(self.k)):
            self.k[i] = self.m.addVar(name="k_"+str(i), vtype = GRB.BINARY)

    def addConstraints(self):
        cnt = 0
        for key in self.angle_diffs_dic:
            self.m.addConstr( self.diff[cnt] == self.a[key[0]] + self.k[cnt] * 360 - self.a[key[1]] - self.angle_diffs_dic[key])
            cnt += 1
        for i in range(len(self.angle_diffs_dic)):
            self.m.addConstr(self.abs_diff[i] == gp.abs_(self.diff[i]))
        for key in self.angle_diffs_constr_dic:
            self.m.addConstr( self.a[key[0]] + self.k[cnt] * 360 - self.a[key[1]], GRB.EQUAL, self.angle_diffs_constr_dic[key])
            cnt += 1

    def solve(self):
        objective = gp.LinExpr()
        for i in range(len(self.angle_diffs_dic)):
            objective += self.abs_diff[i]

        # Set objective: maximize x
        self.m.setObjective(objective, GRB.MINIMIZE)

        try:
            self.m.write("m.lp")
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

opt = Opt(angle_num, angle_diffs_constr_dic, angle_diffs_dic)
opt.createVars()
opt.addConstraints()
opt.solve()
