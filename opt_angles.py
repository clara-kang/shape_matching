import gurobipy as gp
from gurobipy import GRB

# angles = [90, 180, 270, 0, 120, 210, 300, 30]
# Create a new model
m = gp.Model("model")

# Create variables
a = [None] * 8
k = [None] * 8
diff = [None] * 8
abs_diff = [None] * 8
kp = None

for i in range(8):
    a[i] = m.addVar(name="a_"+str(i), ub = 360)
    k[i] = m.addVar(name="k_"+str(i), vtype = GRB.BINARY)
    diff[i] = m.addVar(name="diff_"+str(i), lb = -360, ub = 360)
    abs_diff[i] = m.addVar(name="abs_diff_"+str(i), ub = 360)
kp = m.addVar(name="kp", vtype = GRB.BINARY)

m.addConstr(a[5] + kp * 360 - a[3], GRB.EQUAL, 180)

for face_id in range(2):
    offset = face_id * 4
    for i in range(4):
        m.addConstr( diff[i + offset] == a[(i+1)%4 + offset] + k[i + offset] * 360 - a[i + offset] - 90 )

for i in range(8):
    m.addConstr(abs_diff[i] == gp.abs_(diff[i]))

objective = gp.LinExpr()
for i in range(8):
    objective += abs_diff[i]

# Set objective: maximize x
m.setObjective(objective, GRB.MINIMIZE)

# First optimize() call will fail - need to set NonConvex to 2
try:
    m.write("m.lp")
    # m.feasRelaxS(0, True, False, True)
    # m.write("feasopt1.lp")
    m.optimize()
except gp.GurobiError:
    print("Optimize failed due to non-convexity")


# f = open("square_pts.txt", "w")
# for i in range(24):
#     print("x_", i, ": ", px[i].x)
#     print("y_", i, ": ", py[i].x)
#     f.write(str(px[i].x) + " " + str(py[i].x) + "\n")
# f.close()
#
for v in m.getVars():
    print('%s %g' % (v.varName, v.x))
