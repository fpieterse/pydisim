'''
This file contains descriptions of processes that can be used to create large
systems by connecting the processes. There is an AbstractProcess from which all
processes derive. The requirement for a Process is to have a _connectionInfo
array that the MainProcess will use to record which downstream processes to
update after executing a processs. A Process must also define a run_for function
that performs calculations when the process is run.

The MainProcess class allows you to connect multiple processes together. The way
this class connects the outputs from one class to another expects the values it
connects to be scalars. This is a limitation because we can't specify an
arbitrary number of input streams for things like the MixerProcess. It will be
nice to define a standard stream class that contains composition, flowrate,
temperature, and pressure but this puts more requirements on the subprocesses to
handle all these properties. The idea is to create simple processes. The
downside is that all processes are not necesarilly compatible.
'''

import numpy # for noise

class AbstractProcess:
    '''
    Abstract process unit.

    Represents a process unit with inputs and outputs that can be simulated for
    a given number of seconds by the MainProcess process unit (that is itself
    also a process unit).

    Inherited classes must implement the run_for-function.
    Inherited classes must execute the super-class' __init__-function that
    creates a list of connection strings used by the MainProcess.

    '''

    def __init__(self):
        # Record of connection info used by MainProcess
        # This will be a list of tuples of the form
        # (outputName, downstreamProcess, downstreamInputName)
        self._connectionInfo = []

    def run_for(self,dt):
        '''
        Run process for dt seconds.
        '''
        # Must be implemented
        raise NotImplementedError

class MainProcess(AbstractProcess):
    '''
    The main process represents a collection of subprocesses and their
    connections and relationships
    '''
    def __init__(self):
        super().__init__()
        self._subProcesses = []

    def add_process(self,process):
        '''
        Adds a process to the list of processes
        '''
        self._subProcesses.append(process)

    def add_connection(self,upstreamProcess,
                           upstreamProcessOutput,
                           downstreamProcess,
                           downstreamProcessInput,
                           init=None):
        '''
        Connects output parameter of one process to another one

        Parameters
        ----------
        upstreamProcess : AbstractProcess
            A subprocess to read results from
        upstreamProcessOutput : string
            Name of output of upstream process
        downstreamProcess : AbstractProcess
            A process to write results to
        downstreamProcessOutput : string
            Name of input of downstream process
        init : string
            A string indicating how to initialise values
            <    : Initialise upstream value to equal downstream value.
                   This is most common used initialisation. Use this to
                   initialise PID controller outputs for example.
            >    : Initialise downstream value to equal upstream value
                   This is not very common because if processes are executed in
                   a reasonable order then the values will be overwritten during
                   the first execution. This can be useful to initialise
                   processes that are part of a recycle.
            None : (Default) no initialisation
        '''
        upstreamProcess._connectionInfo.append(
            (upstreamProcessOutput,
             downstreamProcess,
             downstreamProcessInput))

        if init == "<":
            setattr(upstreamProcess,upstreamProcessOutput,
                    getattr(downstreamProcess,downstreamProcessInput))
        elif init == ">":
            setattr(downstreamProcess,downstreamProcessInput,
                    getattr(upstreamProcess,upstreamProcessOutput))


    def run_for(self,dt):
        '''
        Run main process for dt seconds.
        This involves running each subprocess (in the order they were added to
        the main process) and then writing that subprocess' outputs to the
        relevant downstream processes.
        '''
        for proc in self._subProcesses:
            proc.run_for(dt)
            for c in proc._connectionInfo:
                setattr(c[1],c[2],getattr(proc,c[0]))


class HoldupProcess(AbstractProcess):
    '''
    Simulated tank.
    Tank level can go beyon 0-100%


    Parameters:
    -----------
    vol : float
        volume of vessel (m3)

    Simulated Inputs:
    -----------------
    fIn : float
        flow into tank (m3/hr)
    fOut : float
        flow out of tank (m3/hr)
        
    Simulated States:
    -----------------
    cVol : float
        volume of tank contents (vol*level/100) (m3)
    level : float
        liquid level in tank (%)

    Simulated Outputs:
    ------------------
    '''

    @property
    def level(self):
        '''
        The level state is internally maintained as cVol. This property allows
        you to set the cVol as a level %
        '''
        return 100.0 * self.cVol/self.vol
    @level.setter
    def level(self,value):
        self.cVol = self.vol * (value/100.0)
        
    def __init__(self):
        super().__init__()

        self.vol = 1.0
        self.cVol = 0.5

        self.fIn = 0
        self.fOut = 0

    def run_for(self,dt):
        self.cVol += dt*(self.fIn - self.fOut)/3600.0

class SepProcess(AbstractProcess):
    '''
    Binary Separator Process

    Simulates separation process. The simulated process contains one input
    stream and two output streams.

    Simulation Parameters:
    ---------------------
    vol      : Vessel volume (m3)
    relVol   : Relative volatility

    Simulation Inputs:
    -----------------
    F_in     : Feed flowrate (m3/h)
    F_bot    : Bottom flowrate (m3/h)
    F_top    : Top flowrate (m3/h)
    xA_in    : Feed fraction of component A

    Simulation States:
    -----------------
    level    : Tank level (%)
    cVol     : Volume of contents in tank
    xA       : Fraction of component A

    Simulation Outputs:
    ------------------
    xA_top   : Component A fraction in top product
    xA_bot   : Component A fraction in bottom product
    '''

    @property
    def level(self):
        return 100.0 * self.cVol/self.vol
    @level.setter
    def level(self,value):
        self.cVol = (value/100)*self.vol
    @property
    def xA_bot(self):
        return self.xA



    def __init__(self):
        super().__init__()
        self.vol = 1.0
        self.relVol = 2.0
        self.Cp = 4200*1000
        
        self.cVol = 0.5
        self.xA = 0.5
        self.xA_top = 0.5

        self.F_in = 0.0
        self.xA_in = 0.5
        self.F_bot = 0.0
        self.F_top = 0.0

    def run_for(self,dt):
        dtHr = dt/3600

        # Total component A in vessel
        Atot = (self.F_in*self.xA_in*dtHr + self.xA*self.cVol)
        
        self.cVol += self.F_in*dtHr
        if (self.cVol > 0):
            xA = Atot/self.cVol
        else:
            # if the vessel is empty and nothing came in then the concentration
            # stays the same
            xA = self.xA

        # Relative volatility:
        # relVol = (ya/xa)/(yb/xb)
        # yb = 1-ya  ;   xb = 1-xa
        # ya = xa*relVol / ( 1 - xa*(1-relVol) )

        # yA if top flow tends to zero
        yA_NoVap = ( (xA * self.relVol)
                    /(1 - xA*(1-self.relVol)) )
        xA_NoVap = xA
        yA_AllVap = xA
        xA_AllVap = xA/(self.relVol*(1-xA) + xA)

        # What fraction of cVol is top flow
        v_top = self.F_top*dtHr
        if (v_top >= self.cVol):
            v_top = self.cVol
            self.F_top = self.cVol/dtHr
            self.F_bot = 0
            self.cVol = 0
            self.xA = xA_AllVap
            self.xA_top = yA_AllVap
        else:
            fx = v_top/self.cVol
            self.xA_top = fx*yA_AllVap + (1-fx)*yA_NoVap

            Atot -= self.xA_top*v_top
            self.cVol -= v_top
            self.xA = max(Atot/self.cVol, 0)
           
            v_bot = self.F_bot*dtHr
            if v_bot >= self.cVol:
                v_bot = self.cVol
                self.F_bot = v_bot/dtHr
                self.cVol = 0
            else:
                self.cVol -= v_bot

class MixerProcess(AbstractProcess):
    '''
    Mixes streams. Calculates the composition and flowrate of mixture.

    Simulation Inputs:
    ------------------
    F1_F    : Stream 1 flowrate
    F1_xA   : Stream 1 fraction component A
    F2_F    : Stream 2 flowrate
    F2_xA   : Stream 2 fraction component A
    F3_F    : Stream 3 flowrate
    F3_xA   : Stream 3 fraction component A

    Simulation Outputs:
    -------------------
    Fout_F  : Output flowrate
    Fout_xA : Output fraction component A
    '''

    def __init__(self):
        super().__init__()

        # It would be nice to do this as an array but the way MainProcess gets
        # the values of subprocessses to pass around expects floats.
        
        self.F1_F = 0
        self.F1_xA = 0.5
        self.F2_F = 0
        self.F2_xA = 0.5
        self.F3_F = 0
        self.F3_xA = 0.5

        self.Fout_F = 0
        self.Fout_xA = 0.5

    def run_for(self,dt):
        self.Fout_F = self.F1_F + self.F2_F + self.F3_F
        self.Fout_xA = ( (  (self.F1_F*self.F1_xA)
                           +(self.F2_F*self.F2_xA)
                           +(self.F3_F*self.F3_xA))
                        /(self.Fout_F) )

    

class PIDProcess(AbstractProcess):
    '''
    PID controller

    Simulated Parameters:
    ---------------------
    K : float
        Proportional Action. Default action is reverse acting (increase in PV
        causes decrease in OP, like a flow controller), set K negative for
        direct acting (e.g. like a level controller writing to the output of a
        tank) (default 1.0)
    Ti : float
        Integral Time (seconds) (default (1000.0)
    Td: float
        Derivative Time (seconds) (default (0.0)
    PVRange : float
        PV range used for scaling (default 100)
    OPRange : float
        OP range used for scaling (default 100)
    OPLimits : (float,float)
        OP output limits  (default (0,100))
    K_onErr : bool
        Calculate proportional action on Error. Set to false to calculate on PV
        (default = false)
    D_onErr : bool
        Calculate derivative action on Error: Set to false to calculate on PV
        (default = false)
    man : bool
        Manual mode (default = false). In manual mode the execution is skipped

    Simulated Inputs:
    -----------------
    sp : float
        Setpoint [EU]
    pv : float
        Process value [EU]
    K : float
        Proportional Action
    Ti : float
        Integral Time
    Td : float
        Derivative Time

    Simulated States:
    -----------------
    op : float
        Output [EU].
        
    Simulated Outputs: 
    ------------------

    '''

    def get_pv(self):
        return self._nextPv
    def set_pv(self,value,init=False):
        '''
        Used by pv property to set the pv. You can call it directly if you
        wish to initialise the error.

        Parameters:
        -----------
        value : float
            new pv
        init : bool
            set to True to make the change in PV zero
        '''
        self._nextPv = value
        if init:
            self._lastPv = value
            self._lastDxDt = 0.0
    pv = property(get_pv,set_pv)

    def get_sp(self):
        return self._nextSp
    def set_sp(self,value,init=False):
        '''
        Used by sp property to set the setpoint. You can call it directly if you
        wish to initialise the error.

        Parameters:
        -----------
        value : float
            new setpoint
        init : bool
            set to True to make the change in setpoint zero
        '''
        self._nextSp = value
        if init:
            self._lastSp = value
            self._lastDxDt = value
    sp = property(get_sp,set_sp)
            

    def __init__(self):
        super().__init__()

        self.K = 1.0
        self.Ti = 1000.0
        self.Td = 0.0
        self.pvRange = 100.0
        self.opRange = 100.0
        self.opLimits = (0.0,100.0)

        self.K_onErr = False
        self.D_onErr = False

        self._firstRun = True

        self._lastPv = 0.0
        self._lastSp = 0.0
        self._nextPv = 0.0
        self._nextSp = 0.0
        self._lastDxDt = 0.0 # rate of change of derivative term

        self.op = 0.0

        self.man = False


    def run_for(self,dt):
        if (self._firstRun):
            self._lastPv = self._nextPv
            self._lastSp = self._nextSp
        dPv = self._nextPv - self._lastPv

        # change in erorr is change in PV minus change in SP
        dErr = dPv - (self._nextSp - self._lastSp)

        self._lastPv = self._nextPv
        self._lastSp = self._nextSp

        if not self.man:
            # For derivative action we must calculate the rate of change of the
            # derivative term (pv or error) and then calculate the change in the
            # rate of change
            if (self.D_onErr):
                dXdt = dErr/dt
            else:
                dXdt = dPv/dt
            dDD = dXdt - self._lastDxDt
            self._lastDxDt = dXdt

            # dOP is scaled delta OP. Reverse acting is the norm, that is why
            # everything is negative

            # Integral Action
            dOP = -dt*(self.K/self.Ti)*(self._nextPv - self._nextSp)/self.pvRange

            # Proportional Action
            if self.K_onErr:
                dOP -= self.K*(dErr/self.pvRange)
            else:
                dOP -= self.K*(dPv/self.pvRange)

            # Derivative Action
            dOP -= self.K*self.Td*dDD/self.pvRange

            self.op += dOP*self.opRange

        self.op = max(self.opLimits[0],min(self.opLimits[1],self.op))
        self._firstRun = False

class OperatorProcess(AbstractProcess):
    '''
    Simulates operator control actions.

    Simulation Parameters:
    ----------------------
    K            : Integral gain
    reactionTime : (min,max) Operator reacts not faster than minimum time and
                   is forced to make a move after maximum time of no
                   controlRange action was required for period.
                   (hours)
    controlRange : (min,max) Operator reacts when pv out of range
    opRange      : (min,max) Control limits

    Simulation States:
    ------------------
    tMove        : time since last action

    Simulation Inputs:
    ------------------
    pv           : Process value

    Simulation Outputs:
    -------------------
    op           : Output

    '''

    def __init__(self):
        super().__init__()
        self.K = 1.0
        self.reactionTime = (0.1, 2.0)
        self.controlRange = (0.0, 1.0)
        self.pv = 0.5
        self.op = 0.5
        self.tMove =0.0
        self.opRange = (0.0, 1.0)

    def run_for(self,dt):
        dtHr = dt/3600

        self.tMove += dtHr

        if self.tMove < self.reactionTime[0]:
            return

        if   ( self.pv < self.controlRange[0] ):
            self.op += self.K*(self.controlRange[0] - self.pv)
            self.tMove = 0.0
        elif ( self.pv > self.controlRange[1] ):
            self.op += self.K*(self.controlRange[1] - self.pv)
            self.tMove = 0.0
        elif ( self.tMove > (numpy.random.uniform()*self.reactionTime[1]) ):
            sp = (self.controlRange[0] + self.controlRange[1])/2
            self.op += 0.5*self.K*(sp - self.pv)
            self.tMove = 0.0

        self.op = max(self.opRange[0],min(self.opRange[1],self.op))

class SplitRangeProcess(AbstractProcess):
    '''
    Split an iput according to a fraction

    Simulation Parameters:
    ----------------------
    R  : Split ratio (0 - 1)

    Simulation Inputs:
    ------------------
    input : Input signal

    Simulation Outputs:
    -------------------
    op_1 : Output signal 1 (R*input)
    op_2 : Output signal 2 ((1-R)*input)
    '''

    @property
    def op_1(self):
        return self.R * self.input
    @op_1.setter
    def op_1(self,value):
        if self.R > 0.0:
            self.input = value/self.R

    @property
    def op_2(self):
        return (1-self.R)*self.input
    @op_2.setter
    def op_2(self,value):
        if self.R < 1.0:
            self.input = value/(1-self.R)


    def __init__(self):
        super().__init__()
        
        self.R = 0.5
        self.input = 0

    def run_for(self,dt):
        pass


        
class GaussNoiseProcess(AbstractProcess):
    '''
    Adds noise to an input to simulate measurement error

    Simulation Parameters:
    ----------------------
    sigma  : Standard deviation of noise

    Simulation Inputs:
    ------------------
    input  : Input

    Simulation Outputs:
    -------------------
    output : Output (noise + input)

    '''

    def __init__(self):
        super().__init__()

        self.sigma = 0.05
        self.input = 0.0
        self.output = 0.0

    def run_for(self,dt):
        self.output = self.input + numpy.random.normal(0,self.sigma)

class BrownNoiseProcess(AbstractProcess):
    '''
    Simulates process noise. Process noise is Brown noise (I think) (I think
    process noise is brown noise and I think what I did here is brown noise, but
    I won't bet money on it.

    Simulation Parameters:
    ----------------------
    limits  : (low, high) Tuple with low and high limits
    rate    : Rate of movement in standard deviation per hour

    Simulation States:
    ------------------
    output  : Output (noise + input)

    '''

    def __init__(self):
        super().__init__()

        self.output = 0.5
        self.rate = 1.0
        self.limits = (0,1.0)

    def run_for(self,dt):
        self.output += numpy.random.normal(0, self.rate*(dt/3600))
        self.output = min(self.limits[1],max(self.limits[0],self.output))





