#!/usr/bin/env python

import sys
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy

import pydisim as pds

opr = pds.OperatorProcess()

opr.run_for(24*3600)
assert opr.op == 0.5
opr.pv = 1.1
opr.tMove = 0.0
opr.run_for(1)
assert opr.op == 0.5
opr.run_for(0.1*3600)
assert abs(opr.op - 0.4) < 0.001

opr.pv = 0.6
opr.tMove = 0.0
opr.op = 0.5
t = []
rec = pds.Recorder(['op'])
ti = 0
dt = 60
for i in range(1000):
    ti += dt
    t.append(ti)
    opr.run_for(dt)
    rec.record(opr.op)

rec.plot()
plt.show()

sr = pds.SplitRangeProcess()
sr.input = 1
sr.R = 0.4
assert sr.op_1 == 0.4
assert sr.op_2 == 0.6
sr.op_1 = 0.8
assert sr.op_1 == 0.8
assert sr.input == 2.0
sr.op_2 = 0.6
assert sr.input == 1.0
sr.R = 0
sr.op_1 = 1000
assert sr.input == 1.0
sr.R = 1.0
sr.op_2 = 1000
assert sr.input == 1.0
