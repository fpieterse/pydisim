#!/usr/bin/env python
'''
Tests for processes
'''

import sys
sys.path.append('../')

import matplotlib.pyplot as plt

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
pid.RunFor(1)
assert pid.op == 50

# Test proportional-only control
pid.Ti = float('inf')
pid.pv = 51
pid.RunFor(1)
assert pid.op == 49

# Test integral action
pid.Ti = 100.0
pid.RunFor(100)
assert pid.op == 48

# Test K_onPv - the setpoint changes so integral action will be zero.
pid.sp = 51
pid.RunFor(100)
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
pid.RunFor(100)
assert pid.op == 2.0 # output increase due to sp-pv error

pid.sp = 51
pid.RunFor(100)
assert pid.op == 0.0 # output decrease due to sp change (no integral action)



#---- TANK -----------------------------------------------------------
tank = HoldupProcess(1,50)
tank.fIn = 1.0
tank.fOut = 1.0
tank.RunFor(3600)
assert tank.level == 50

tank.fOut = 1.1
tank.RunFor(3600)
assert abs(tank.level - 40) < 0.0000001


#---- Tank level control -----------------------------------------------------
# Run a tank with a level controller with fixed tuning for a while and confirm
# that the final state is at setpoint. I also track the turning points and
# maximum level and flow output to confirm this works as expected
#-----------------------------------------------------------------------------
tk001 = HoldupProcess(1,50)
tk001.fIn = 1.5
tk001.fOut = 1.0

lc001 = PIDProcess()
lc001.K = -1
lc001.Ti = 3600
lc001.opRange = 2.0
lc001.opLimits = (0,2.0)
lc001.sp = 50

proc = MainProcess()
proc.AddProcess( tk001 )
proc.AddProcess( lc001 )
proc.AddConnection(tk001,'level',lc001,'pv')
proc.AddConnection(tk001,'fOut',lc001,'op')
proc.AddConnection(lc001,'op',tk001,'fOut')

t = 0
dt = 60
maxOp = 0
maxOpT = 0
maxLvl = 0.0
maxLvlT = 0
for i in range(500):
    t += dt
    proc.RunFor(dt)

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


#---- Reactors -----------------------------------------------------------
# CSTR with no inflow
rx = CSTRProcessI()
rx.T = 230
rx.FJ_T = 300
rx.FJ_F = 1.0

time = []
T = []
TJ = []
xA = []
xB = []
xC = []

dt = 10
t = 0
for i in range(10):
    t += dt
    time.append(t)
    rx.RunFor(dt)
    T.append(rx.T)
    TJ.append(rx.J_T)
    xA.append(rx.xA)
    xB.append(rx.xB)
    xC.append(rx.xC)

plt.plot(time,xA)
plt.plot(time,xB)
plt.plot(time,xC)
plt.show()

plt.plot(time,T)
plt.plot(time,TJ)
plt.show()


