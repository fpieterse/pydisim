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

main = pds.MainProcess()
bnp = pds.BrownNoiseProcess()
gnp = pds.GaussNoiseProcess()
fp = pds.FilterProcess()
fp.t = f

main.add_process(bnp)
main.add_process(gnp)
main.add_process(fp)

main.add_connection(bnp,'output',gnp,'input')
main.add_connection(gnp,'output',fp,'input')

rec = pds.tools.Recorder(['Brown','Gauss','Filter'])

# PLot random noise. Then do it again (need to see different noise)
# THen do it again with seed to make sure we can seed the noise with numpy

print("This is noise, the next plot must look different than this one")
for i in range(100):
    main.run_for(10)
    rec.record(bnp.output,gnp.output,fp.output)

rec.plot()
plt.show()

rec.clear()

print("This is noise, it must be different than the first one but same as the next")
numpy.random.seed(0)
for i in range(100):
    main.run_for(10)
    rec.record(bnp.output,gnp.output,fp.output)

rec.plot()
plt.show()
rec.clear()

print("This noise must look the same as the previous one")
numpy.random.seed(0)
for i in range(100):
    main.run_for(10)
    rec.record(bnp.output,gnp.output,fp.output)

rec.plot()
plt.show()

