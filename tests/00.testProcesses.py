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

#---- Deadtime ---------------------------------------------------

dt = DeadtimeProcess(deadtime=10)
dt.output = 5
assert dt.input == 5

dt.input = 3
assert dt.history == [3]
assert dt.t == [0]

dt.run_for(1)
assert dt.t == [1,0]

dt.input = 2
dt.run_for(2)
assert dt.t == [3,2,0]


dt.run_for(7)
assert dt.t == [10,9,7,0]
assert dt.history == [3,2,2,2]
assert dt.output == 3

dt.run_for(1)
dt.run_for(1)
dt.run_for(1)
assert dt.output == 2

dt.input = 5
dt.run_for(1)
assert dt.output == 2


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


# Test OP does not wind up
pid.set_pv(51,True)
pid.set_sp(50,True)
pid.op = 100
pid.pvRange = 100
pid.opRange = 100
pid.opLimits = (0,100)
pid.K = -1.0
pid.Ti = 1
pid.run_for(1)
assert pid.co == 100
pid.pv += 1
pid.run_for(1)
assert pid.co == 101
assert pid.op == 100
pid.pv -= 1
pid.run_for(1)
assert pid.co == 100




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

sep.F_in = 1.0
sep.level = 50
sep.F_top = 0.5
sep.F_bot = 0.5
sep.run_for(1000)
assert abs(sep.level - 50) < 0.001
sep.F_top = 0.4
sep.run_for(3600)
assert abs(sep.level -60) < 0.001
sep.F_bot = 0.6
sep.run_for(3600)
assert abs(sep.level -60) < 0.001

#---- MIXER-------------------------------------------------------------------

mix = MixerProcess()
mix.F1_F = 1.0
mix.F2_F = 3.0
mix.F3_F = 4.0
mix.F1_xA = 0.0
mix.F2_xA = 1.0
mix.F3_xA = 0.5
mix.run_for(1)
assert mix.Fout_F == 8
assert mix.Fout_xA == 0.625


#---- SELECTOR----------------------------------------------------------------

sel = SelectProcess(n_inputs=5)
assert len(sel.input) == 5
sel.input[0] = 1
sel.input[1] = 99
sel.input[2] = -10
sel.input[3] = 8
sel.input[4] = 5

sel.seltype = 'high'
sel.run_for(1)
assert sel.output == 99

sel.seltype = 'low'
sel.run_for(1)
assert sel.output == -10

sel.seltype = 'median'
sel.run_for(1)
assert sel.output == 5

# setting sel.output will overwrite the inputs
sel._output = 0
sel.rate = 1
sel.run_for(1)
assert sel.output == 1


#---- MATHADD MATHMUL --------------------------------------------------------
add = MathAddProcess(n_inputs=5)
add.bias = 11.5
add.scale[0] = 1
add.scale[1] = 2
add.scale[2] = 44.5
add.scale[3] = 0
add.scale[4] = -10

add.output = 10
assert add.output == 10

# if all the scales are zero, the inputs should be zero
add.scale[:] = 0
add.output = 10
assert add.output == add.bias 
assert add.input[0] == 0

# if the sum of the scales is zero then the inputs cannot be calculated
# confirm no errors are raised and the output is equal to the bias.
add.scale[1] = 1
add.output = 10
assert add.output == 10
add.scale[2] = -1
add.output = 10
assert add.output == 10


mul = MathMulProcess(n_inputs=3)
mul.bias = 1
mul.scale = 10
mul.output = 10
assert mul.output == 10
mul.scale = 0
mul.output = 10
assert mul.output == 10



#---- Tank level control -----------------------------------------------------
# Run a tank with a level controller with fixed tuning for a while and confirm
# that the final state is at setpoint. I also track the turning points and
# maximum level and flow output to confirm this works as expected
#-----------------------------------------------------------------------------

# Create a new process manager
pm = ProcessManager()

tk001 = HoldupProcess()
tk001.fIn = 1.5
tk001.fOut = 0.0 # will overwrite in initialisation

tk002 = HoldupProcess()
tk002.fIn = 1.0
tk002.fOut = 1.0

lc001 = PIDProcess()
lc001.K = -1
lc001.Ti = 3600
lc001.opRange = 2.0
lc001.opLimits = (0,2.0)
lc001.sp = 50

lc002 = PIDProcess()
lc002.K = -1
lc002.Ti = 1800
lc002.opRange = 2.0
lc002.opLimits = (0,2.0)
lc002.sp = 50

tk001.add_connection('fOut',tk002,'fIn','<')
tk001.add_connection('level',lc001,'pv','>')
lc001.add_connection('op',tk001,'fOut','<')
tk002.add_connection('level',lc002,'pv')
lc002.add_connection('op',tk002,'fOut','<')

assert tk001.fOut == 1.0
assert lc001.pv == 50.0
assert lc001.op == 1.0
assert lc002.pv == 0.0
assert lc002.op == 1.0

t = 0
# First test ran process at 60 second steps and checked at which point in the
# execution the maximum value is observed (if dynamics are not affected then the
# point should be the same).  With the introduction of the ProcessManager, I now
# execute the process manager at 60 sec scan time and run it for one minute at a
# time (so it does one execution).
pm.exec_int = 60
maxOp = 0
maxOpT = 0
maxLvl = 0.0
maxLvlT = 0
for i in range(500):
    t += dt
    #Run process for one minute (at default 5 sec interval)
    pm.run_process(1)

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


