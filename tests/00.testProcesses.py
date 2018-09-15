#!/usr/bin/env python
'''
Tests for processes
'''

import sys
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy
import scipy.interpolate

from pydisim.processes import *

#---- PID CONTROLLER ---------------------------------------------

pid =  PIDProcess()
assert pid.K == 1.0
assert pid.Ti == 1000.0
assert pid.Td == 0.0

# The first time we run the PID it should do nothing
pid.op = 50
pid.pv = 50
pid.sp = 50
pid.run_for(1)
assert pid.op == 50

# Test proportional-only control
pid.Ti = float('inf')
pid.pv = 51
pid.run_for(1)
assert pid.op == 49

# Test integral action
pid.Ti = 100.0
pid.run_for(100)
assert pid.op == 48

# Test K_onPv - the setpoint changes so integral action will be zero.
pid.sp = 51
pid.run_for(100)
assert pid.op == 48

# Test scaling and other advanced functions
pid.set_pv(51,True)
pid.set_sp(50,True)
pid.pvRange = 10.0
pid.opRange = 20.0
pid.op = 0
pid.opLimits = (-10,10)
pid.K_onErr = True
pid.K = -1.0
pid.Ti = 100.0
pid.run_for(100)
assert pid.op == 2.0 # output increase due to sp-pv error

pid.sp = 51
pid.run_for(100)
assert pid.op == 0.0 # output decrease due to sp change (no integral action)



#---- TANK -----------------------------------------------------------
tank = HoldupProcess()

assert tank.level == 50
tank.level = 60
assert tank.cVol == 0.6
tank.fIn = 1.0
tank.fOut = 1.0
tank.run_for(3600)
assert tank.level == 60

tank.fOut = 1.1
tank.run_for(3600)
assert abs(tank.level - 50) < 0.0000001


#---- SEP --------------------------------------------------------------------
sep = SepProcess()

assert sep.level == 50

# Liquid Vapour Equilibrium
sep.relVol = 3.0
x_vle = numpy.linspace(0,1)
y_vle = []
for i in range(len(x_vle)):
    sep.xA = x_vle[i]
    sep.run_for(1)
    y_vle.append( sep.xA_top )

y_vle = numpy.array(y_vle)

# Batch distillation test
sep.xA = 0.9
sep.F_top = 0.1

yA = []
xA = []
t = []
ti = 0
dt = 60
while (sep.F_top > 0.01):
    ti += dt
    t.append(ti)
    sep.run_for(dt)
    xA.append(sep.xA)
    yA.append(sep.xA_top)

#plt.plot(xA,yA)
#plt.plot(x_vle,y_vle)
#plt.show()
vleInterp = scipy.interpolate.interp1d(x_vle,y_vle)
y_check = vleInterp(xA)
yA = numpy.array(yA)
assert (sum(abs(y_check - yA))/len(yA)) < 0.01
print(sum(abs(y_check - yA))/len(yA))

# Full vap / full liq test
sep.xA = 0.5
sep.level = 50
assert sep.cVol == 0.5

sep.F_top = 0.0
sep.run_for(1)
y_NoVap = sep.xA_top
assert sep.xA_top == 0.75
sep.F_top = 0.5
sep.run_for(3600)
y_AllVap = sep.xA_top
assert sep.xA_top == 0.5

# Take 50% of volume out, confirm the top product is purer
# but less pure than when there is no vapour
sep.xA = 0.5
sep.level = 50
sep.F_top = 0.25
sep.run_for(3600)
assert sep.xA_top == 0.625
assert sep.xA_bot == 0.375


#---- Tank level control -----------------------------------------------------
# Run a tank with a level controller with fixed tuning for a while and confirm
# that the final state is at setpoint. I also track the turning points and
# maximum level and flow output to confirm this works as expected
#-----------------------------------------------------------------------------
tk001 = HoldupProcess()
tk001.fIn = 1.5
tk001.fOut = 1.0

lc001 = PIDProcess()
lc001.K = -1
lc001.Ti = 3600
lc001.opRange = 2.0
lc001.opLimits = (0,2.0)
lc001.sp = 50

proc = MainProcess()
proc.add_process( tk001 )
proc.add_process( lc001 )
proc.add_connection(tk001,'level',lc001,'pv')
proc.add_connection(tk001,'fOut',lc001,'op')
proc.add_connection(lc001,'op',tk001,'fOut')

t = 0
dt = 60
maxOp = 0
maxOpT = 0
maxLvl = 0.0
maxLvlT = 0
for i in range(500):
    t += dt
    proc.run_for(dt)

    if (tk001.level > maxLvl):
        maxLvlT = t
        maxLvl = tk001.level
    if (lc001.op > maxOp):
        maxOpT = t
        maxOp = lc001.op

#print('MaxOpT:',maxOpT,'MaxOP',maxOp,'MaxLvl',maxLvl,'maxLvlT',maxLvlT)
assert maxOpT == 5580
assert maxLvlT == 2820
assert abs(1.6 - maxOp) < 0.1
assert abs(66.7 - maxLvl) < 0.1
assert abs(50 - tk001.level) < 0.1


