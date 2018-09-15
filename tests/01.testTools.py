#!/usr/bin/env python
'''
Tests for tools
'''

class TestPID():
    '''
    This class is just used as a data container to test the pid recorder
    '''
    sp = 0
    pv = 0
    op = 0

doPlot = False

import sys
sys.path.append('../')

if doPlot:
    import matplotlib.pyplot as plt

from pydisim.tools import *


t = []
var1 = 1.0

rec = Recorder(['Variable1','Variable2'])

for i in range(9):
    t.append(i)
    rec.record(i,i*2,3)

assert rec.data[0][0] == 0
assert rec.data[1][8] == 16

rec.clear(3)
assert len(rec.data) == 2
assert len(rec.data[0]) == 3
assert rec.data[0][2] == 8


if doPlot:
    rec.plot(plt,t)
    plt.show()

rec = Recorder(['v1','v2'])
prec = PIDRecorder('p')
rec.record(1,10)
rec.record(2,20)
pid = TestPID()
pid.sp = 3
pid.pv = 4
pid.op = 5
prec.record(pid)
pid.sp = 6
pid.pv = 7
pid.op = 8
prec.record(pid)
time = [1,2]
savefile = '01.testTools.savetest.csv'
save_to_file(savefile,time,rec,prec)
with open(savefile,'r') as f:
    assert f.readline() == 'time,v1,v2,p.sp,p.pv,p.op\n'
    assert f.readline() == '1,1.000000,10.000000,3.000000,4.000000,5.000000\n'
    assert f.readline() == '2,2.000000,20.000000,6.000000,7.000000,8.000000\n'
