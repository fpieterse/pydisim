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

#-----------------------------------------------------
# Sin
#-----------------------------------------------------
pm = pds.ProcessManager(exec_int=1)
swp = pds.SinNoiseProcess(amplitude=5,period=120)
rec = pds.RecorderProcess(['Input','Output'])

swp.add_connection('input',rec,'Input')
swp.add_connection('output',rec,'Output')

# Check that inputs and outputs initiallises correctly
swp.output = 5
assert swp.input == 5
swp.input = 0
assert swp.output == 0

# check that sin wave works
pm.run_process(0.5) #run for 30 seconds
sin30 = numpy.sin(0.5*numpy.pi)
print(sin30)
print(swp.output)
assert swp.output == sin30
swp.input = 1
assert swp.output == sin30+1

# check that if we change the period the output doesn't bump
swp.period = 60
assert swp.output == sin30+1
pm.run_process(1) #Run for 1 minute
assert swp.output == sin30+1

rec.plot()
plt.show()



#-----------------------------------------------------
# Noise and Filters
#-----------------------------------------------------
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

