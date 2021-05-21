#!/usr/bin/env python

import sys
sys.path.append('../')

f = 20.0
if len(sys.argv) > 1:
    try:
        f = float(sys.argv[1])
    except:
        pass

import matplotlib.pyplot as plt
import numpy

import pydisim as pds

pm = pds.ProcessManager()
bnp = pds.BrownNoiseProcess()
gnp = pds.GaussNoiseProcess()
fp = pds.FilterProcess()
fp.t = f


bnp.add_connection('output',gnp,'input')
gnp.add_connection('output',fp,'input')

rec = pds.tools.RecorderProcess(['Brown','Gauss','Filter'])
bnp.add_connection('output',rec,'Brown')
gnp.add_connection('output',rec,'Gauss')
fp.add_connection('output',rec,'Filter')

# PLot random noise. Then do it again (need to see different noise)
# THen do it again with seed to make sure we can seed the noise with numpy

print("Figure 1 is noise, the next plot must look different than this one")
pm.run_process(5)

rec.plot()
rec.clear()

print("Figure 2 is also noise, it must be different than the first one but same as the next")
numpy.random.seed(0)
pm.run_process(5)

rec.plot()
rec.clear()

print("Figure 3 must look the same as Figure 2.")
numpy.random.seed(0)
pm.run_process(5)

rec.plot()
plt.show()

