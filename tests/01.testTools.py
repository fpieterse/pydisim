#!/usr/bin/env python
'''
Tests for tools
'''

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
prec.record(3,4,5)
prec.record(6,7,8)
time = [1,2]
savefile = '01.testTools.savetest.csv'
save_to_file(savefile,time,rec,prec)
with open(savefile,'r') as f:
    assert f.readline() == 'time,v1,v2,p.sp,p.pv,p.op\n'
    assert f.readline() == '1,1.000000,10.000000,3.000000,4.000000,5.000000\n'
    assert f.readline() == '2,2.000000,20.000000,6.000000,7.000000,8.000000\n'
