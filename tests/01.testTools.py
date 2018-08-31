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

if doPlot:
    rec.plot(plt,t)
    plt.show()
