# pydisim
Discrete process simulation for process control simulation.

# Using pydimis

copy the folder to your working folder.

import pydisim as pds

TODO: Install with PYPI.


# Example
```
import pydisim as pds
import numpy

# STEP 1: Create a Custom process by inheriting AbstractProcess
# Or use a built-in processses from processes.py
class PumpProcess(pds.AbstractProcess):
    '''
    VSD Pump
    Set pump speed 0-100
    '''
    
    def __init__(self):
        super().__init__()
        
        self.max_flow = 10
        self.max_current = 100
        
        self.speed = 0
        self.current = 0
        self.flow = 0
        
        # Dynamic Behaviour (time constant)
        self.tau = 4
        
    def run_for(self,dt):
        if self.tau > 0:
            f = numpy.exp(-dt/self.tau)
        else:
            f = 0.0
            
        if self.speed > 100:
            self.speed = 100
        if self.speed < 0:
            self.speed = 0
         
        # steady-state current
        ss_current = self.speed * self.max_current / 100.0
        # Go To SS
        self.current = self.current*f + ss_current*(1-f)
        
        # stead-state flow
        ss_flow = self.current * self.max_flow / self.max_current
        self.flow = self.flow*f + ss_flow*(1-f)

# STEP 2: create a management object, the management object will execute all the function blocks that are created after
# the management object is created
# function blocks are executed in the order they are created
pm = pds.ProcessManager(exec_int =1)

# STEP 3: Create function blocks (a.k.a processes).
     
# Built-in PID controller
FC201 = pds.PIDProcess()
FC201.K = 1
FC201.Ti = 10
FC201.pvRange = 10
FC201.opRange = 100

FC201.sp = 5

# Instance of custom VSD pump class
pump = PumpProcess()
# You could write Process classes to initialise to a steady state on the first execution
# but I'm just going to guess some values:
pump.speed = 40
pump.flow = 4

# Simulated measurement Noise
meas_current = pds.GaussNoiseProcess()
meas_current.sigma = 3

meas_flow = pds.GaussNoiseProcess()
meas_flow.sigma = 0.5

# PV Filter for flow controller
filt_flow = pds.FirstOrderProcess()
filt_flow.t = 3

# Addition Operator to add a bias to the measured current
bias_current = pds.MathAddProcess(n_inputs=2)


# STEP 4: Connect processes
# When management object executes each process, it will write all the
# connections defined on that process.
# The connection initialisation option is used to syncronise starting values, but is not strictly necessary
# .add_connection( output name,
#                  downstream process,
#                  downstream process input name,
#                  initialisation direction )
FC201.add_connection('op',pump,'speed','<')
pump.add_connection('flow',meas_flow,'input','>')
pump.add_connection('current',meas_current,'input','>')
meas_flow.add_connection('output',filt_flow,'input','>')
filt_flow.add_connection('output',FC201,'pv','>')
meas_current.add_connection('output',bias_current,'input[1]','>')

# A convenient plotting object to plot PIDProcess parameters
# Also see the pds.RecorderProcess for custom plotting
rec = pds.PIDRecorderProcess('FC201',FC201)


# STEP 5: Run the management object
pm.run_process(10)
rec.plot()

# STEP 6..n: Make changes, run again, etc.
FC201.sp = 6
pm.run_process(10)
rec.plot()
```
