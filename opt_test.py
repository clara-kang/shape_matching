#!/usr/bin/env python3.7

# Copyright 2020, Gurobi Optimization, LLC

# This example formulates and solves the following simple bilinear model:
#  maximize    x
#  subject to  x + y + z <= 10
#              x * y <= 2         (bilinear inequality)
#              x * z + y * z = 1  (bilinear equality)
#              x, y, z non-negative (x integral in second version)

import gurobipy as gp
from gurobipy import GRB

# Create a new model
m = gp.Model("bilinear")

# Create variables
px = [None] * 4
py = [None] * 4
vx = [None] * 2
vy = [None] * 2
abs_vx = [None] * 2
abs_vy = [None] * 2

for i in range(4):
    px[i] = m.addVar(name="x_"+str(i))
    py[i] = m.addVar(name="y_"+str(i))

for i in range(2):
    vx[i] = m.addVar(name="vx_"+str(i))
    vy[i] = m.addVar(name="vy_"+str(i))
    abs_vx[i] = m.addVar(name="abs_vx_"+str(i))
    abs_vy[i] = m.addVar(name="abs_vy_"+str(i))

# define vs
m.addConstr(vx[0] == px[1] - px[0])
m.addConstr(vy[0] == py[1] - py[0])
m.addConstr(vx[1] == px[3] - px[2])
m.addConstr(vy[1] == py[3] - py[2])

# v0 v1 orthogonal
m.addConstr(vx[0] * vx[1] + vy[0] * vy[1] == 0)

# define abs x y
for i in range(2):
    m.addConstr(abs_vx[i] == gp.abs_(vx[i]))
    m.addConstr(abs_vy[i] == gp.abs_(vy[i]))

# v length
m.addConstr(abs_vx[0] + abs_vy[0] == 1)
m.addConstr(abs_vx[1] + abs_vy[1] == 1)


# # Set objective: maximize x
# m.setObjective(objective, GRB.MINIMIZE)

# First optimize() call will fail - need to set NonConvex to 2
try:
    m.params.NonConvex = 2
    m.optimize()
except gp.GurobiError:
    print("Optimize failed due to non-convexity")

for v in m.getVars():
    print('%s %g' % (v.varName, v.x))
