#!/usr/bin/env python

import sys
sys.path.append('../')

import matplotlib.pyplot as plt
import numpy

import pydisim as pds

main = pds.MainProcess()
bnp = pds.BrownNoiseProcess()
gnp = pds.GaussNoiseProcess()

main.add_process(bnp)
main.add_process(gnp)

main.add_connection(bnp,'output',gnp,'input')

rec = pds.tools.Recorder(['Brown','Gauss'])

# PLot random noise. Then do it again (need to see different noise)
# THen do it again with seed to make sure we can seed the noise with numpy
for i in range(100):
    main.run_for(10)
    rec.record(bnp.output,gnp.output)

rec.plot()
plt.show()

rec.clear()
numpy.random.seed(0)
for i in range(100):
    main.run_for(10)
    rec.record(bnp.output,gnp.output)

rec.plot()
plt.show()
rec.clear()
numpy.random.seed(0)
for i in range(100):
    main.run_for(10)
    rec.record(bnp.output,gnp.output)

rec.plot()
plt.show()

