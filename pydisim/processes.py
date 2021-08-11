'''
This library component contains descriptions of processes that can be used to
create large systems of connected processes.

The libary defines an AbstractProcess from which all
processes derive. You can create custom processes by subclassing the
AbstractProcess class and defining the run_for function.  Some simple processes
are also defined.

The ProcessManager class will run all the processes in the order they were
created and manage connections between the processes.
'''

import numpy # for noise
import re
import datetime
import scipy


def _write_connection(upstream_process,upstream_name,upstream_index,
                      downstream_process,downstream_name,downstream_index):
    '''
    Write an attribute from upstream process to downstream process.

    This function is used by the AbstractProcess to update downstream processes
    during execution and it is used to initialise connections when connection
    are created (which sometimes go in the other direction).

    Parameters:
    -----------
    upstream_process : AbstractProcess
        Source process of connection

    upstream_name : string
        Name of attribute that is source of connection

    upstream_index : integer, None
        If upstream attribute is a list, upstream_index indicates the index
        of the source attribute to write to the downstream process.
        If upstream attribute is scalar then upstream_index must be None.

    downstream_process: AbstractProcess
        Sink process of connection

    downstream_name : string
        Name of attribute that is sink of connection

    downstream_index : integer, None
        See upstream_index.

    '''

    out_value = getattr(upstream_process,upstream_name)
    if upstream_index != None:
        out_value = out_value[upstream_index]

    if downstream_index != None:
        in_value = getattr(downstream_process,
                           downstream_name)
        in_value[downstream_index] = out_value
    else:
        setattr(downstream_process,
                downstream_name,
                out_value)

_last_process_manager = None
def get_process_manager():
    '''
    Return the last process manager created.
    Creates a new one if it has not been created

    '''
    global _last_process_manager
    if _last_process_manager == None:
        _last_process_manager = ProcessManager()

    return _last_process_manager

class ProcessManager:
    '''
    Class to manage subprocesses.

    '''
    def __init__(self,exec_int=5):
        '''
        Parameters:
        -----------
        exec_int : integer
            Duration of process simulation iterations in seconds.
        '''
        self.exec_int = exec_int
        self._subprocesses = []
        self._recorders = []

        self.simulation_start = datetime.datetime.now()
        self.simulation_time = 0

        global _last_process_manager
        _last_process_manager = self

    def run_process(self,duration):
        '''
        Parameters:
        -----------
        duration : int
            Duration to run process for in minutes
        '''
        t_end = self.simulation_time + 60*duration
        while True:
            self.simulation_time += self.exec_int
            for p in self._subprocesses:
                p._execute(self.exec_int)

            if self.simulation_time >= t_end:
                break


    def clear_recorders(self,keep=None):
        for rec in self._recorders:
            rec.clear(keep)
            

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

    _array_num_re = re.compile('(.*)\[([0-9]*)\]$')

    def __init__(self,process_manager=None):
        # Record of connection info.
        # This is a list of dictionaries that define connection information
        self._connectionInfo = []

        if process_manager == None:
            process_manager = get_process_manager()

        process_manager._subprocesses.append(self)
        self.process_manager = process_manager

    def _execute(self,dt):
        self.run_for(dt)

        # TODO: Write a test for this
        for ci in self._connectionInfo:
            _write_connection(self,ci['output_name'],ci['output_index'],
                              ci['downstream_process'],
                              ci['input_name'], ci['input_index'])

    def run_for(self,dt):
        '''
        Run process for dt seconds.

        Must be implemented on inheriting classes
        
        '''
        raise NotImplementedError


    def add_connection(self,
                       output_name,
                       downstream_process,
                       input_name,
                       init=None):
        '''
        Add a connection to a dowstream process. When this process is run, the
        results will be updated in the downstream process.

        Parameters:
        -----------
        output_name : string
            Name of output to connect, output is accessed using getattr
            function.
            If the output is an array, an element can be specified e.g.
            "output[0]".
            If the input and outputs are lists or tuples they will be set equal.
        downstream_process : AbstractProcess
            Downstream process to connect to
        input_name : string
            Name of parameter on downstream process to write output to.
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
        output_index = None
        input_index = None

        m = self._array_num_re.match(output_name)
        if m:
            output_name = m.group(1)
            output_index = int(m.group(2))
        m = self._array_num_re.match(input_name)
        if m:
            input_name = m.group(1)
            input_index = int(m.group(2))

        _output = getattr(self,output_name)
        _input = getattr(downstream_process,input_name)
        
        # If both indexes are None, check if it arrays, if so set them equal and
        # return.
        # If index not None (a number) and the value is out of range, raise
        # exception

        if (output_index == None) and (input_index == None):
            if ((type(_output) == list) and (type(_input) == list)) or \
               ((type(_output) == tuple) and (type(_input) == tuple)):
                if len(_output) == len(_input):
                    if init == '<':
                        setattr(self,output_name,_input)
                    else:
                        setattr(downstream_process,input_name,_output)
                    return
                else:
                    raise Exception("Cannot connect lists of different lengths")
            elif (type(_output) == numpy.ndarray) and (type(_input) == numpy.ndarray):
                if _output.shape == _input.shape:
                    if init == '<':
                        setattr(self,output_name,_input)
                    else:
                        setattr(downstream_process,input_name,_output)
                    return
                else:
                    raise Exception("Cannot connect ndarrays of different shapes")
        elif (output_index != None) and (output_index >= len(_output)):
            raise Exception("Index out of range on output " + output_name)
        elif (input_index != None) and (input_index >= len(_input)):
            raise Exception("Index out of range on input " + input_name)


        if init == "<":
            _write_connection(downstream_process,input_name,input_index,
                              self,output_name,output_index)
        elif init == ">":
            _write_connection(self,output_name,output_index,
                              downstream_process,input_name,input_index)

        conn_info = {'downstream_process':downstream_process,
                     'output_name':output_name,
                     'input_name':input_name,
                     'output_index':output_index,
                     'input_index':input_index,
                    }

        self._connectionInfo.append(conn_info)

            


class HoldupProcess(AbstractProcess):
    '''
    Simulated tank.
    Tank level can go beyon 0-100%


    Parameters:
    -----------
    vol : float
        volume of vessel (EU)

    Simulated Inputs:
    -----------------
    fIn : float
        flow into tank (EU/hr)
    fOut : float
        flow out of tank (EU/hr)
        
    Simulated States:
    -----------------
    cVol : float
        volume of tank contents (vol*level/100) (EU)
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

        if (self.Fout_F > 0):
            self.Fout_xA = ( (  (self.F1_F*self.F1_xA)
                               +(self.F2_F*self.F2_xA)
                               +(self.F3_F*self.F3_xA))
                            /(self.Fout_F) )
        else:
            self.Fout_xA = 0.5



class FirstOrderProcess(AbstractProcess):
    '''
    First Order Process
    -------------------

    Simulated Parameters:
    ---------------------
    t : Filter time (seconds)

    Simulated Inputs:
    -----------------
    input : float
        Process value [EU]

    Simulated States:
    -----------------
    output : float
        Filtered process value [EU]

    '''

    @property
    def input(self):
        return self._input
    @input.setter
    def input(self,value):
        if (self._firstRun):
            self._output = value
            self._firstRun = False
        
        self._input = value

    @property
    def output(self):
        return self._output
    @output.setter
    def output(self,value):
        if (self._firstRun):
            self._input = value
            self._firstRun = False

        self._output =value

    def __init__(self):
        super().__init__()

        self.t = 0.0
        self._firstRun = True
        
        self._output = 0.0
        self._input = 0.0

    def run_for(self,dt):
        if (self._firstRun):
            self._output = self._input
        self._firstRun = False

        if self.t > 0.0:
            f = numpy.exp(-dt/self.t)
        else:
            f = 0.0
        self._output = f*self._output + (1-f)*self._input

# First order process used to be called FilterProcess
FilterProcess = FirstOrderProcess


class DeadtimeProcess(AbstractProcess):
    '''
    DeadtimeProcess
    ---------------

    Simulated Inputs:
    -----------------
    input : float
    
    Simulated States:
    -----------------
    output : float


    Simulation Parameters:
    ----------------------
    deadtime : float
        Deadtime in seconds
    '''

    @property
    def output(self):
        return self.history[0]
    @output.setter
    def output(self,value):
        self.history = [value]
        self.t = [0]

    @property
    def input(self):
        return self.history[-1]
    @input.setter
    def input(self,value):
        self.history[-1] = value


    def __init__(self,deadtime=0):
        super().__init__()
        self.deadtime = deadtime

        # previous input values
        self.history = [0]
        # timestamps of previous values
        self.t = [0]

    def run_for(self,dt):
        i = 1
        while i <= len(self.t):
            self.t[-i] += dt
            if self.t[-i] >= self.deadtime:
                self.t = self.t[-i:]
                self.history = self.history[-i:]
                break
            i += 1
        self.history.append(self.input)
        self.t.append(0)


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
    pvRange : float
        PV range used for scaling (default 100)
    opRange : float
        OP range used for scaling (default 100)
    opLimits : (float,float)
        OP output limits  (default (0,100))
    oplo : float
        OP Low Limit. Sets the first element in opLimits
    ophi : float
        OP high Limit. Sets the second element in opLimits
    K_onErr : bool
        Calculate proportional action on Error. Set to false to calculate on PV
        (default = false)
    D_onErr : bool
        Calculate derivative action on Error: Set to false to calculate on PV
        (default = false)
    man : bool
        Manual mode (default = false). In manual mode the execution is skipped
    pvtrack : bool
        If true, sp is set equal to pv if man is true.
    non_interacting : bool
        If true, use non-interacting form (Ki = 1/Ti, Kd = 1*Td)
    reverse : bool
        Reverse acting control; when the PV increases, the output decreases.  Set to False for Direct acting. (default = True 
    ff_gain : float
        Feedforward gain

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
    ff : float
        Feedforward Input

    Simulated States:
    -----------------
    op : float
        Output [EU].
        

    '''
    # TODO: Non-Interacting is the wrong word.  Interacting means d-action is
    # multiplied by P+I, the difference I'm using is between the Ideal and
    # Series equations.

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

    @property
    def oplo(self):
        return self.opLimits[0]
    @oplo.setter
    def oplo(self,value):
        self.opLimits = (value, self.opLimits[1])

    @property
    def ophi(self):
        return self.opLimits[1]
    @ophi.setter
    def ophi(self,value):
        self.opLimits = (self.opLimits[0],value)

    @property
    def K(self):
        return self._K
    @K.setter
    def K(self,value):
        if value < 0:
            self._K = -value
            self.reverse = False
        else:
            self._K = value

    @property
    def ff(self):
        return self._nextFF
    @ff.setter
    def ff(self,value,init=False):
        self._nextFF = value
        if init:
            self._lastFF = value
            

    def __init__(self):
        super().__init__()

        self._K = 1.0
        self.Ti = 1000.0
        self.Td = 0.0
        self.pvRange = 100.0
        self.opRange = 100.0
        self.opLimits = (0.0,100.0)
        self.ff_gain = 1

        self.K_onErr = False
        self.D_onErr = False

        self._firstRun = True

        self._lastPv = 0.0
        self._lastSp = 0.0
        self._lastFF = 0.0
        self._nextPv = 0.0
        self._nextSp = 0.0
        self._nextFF = 0.0
        self._lastDxDt = 0.0 # rate of change of derivative term

        self.op = 0.0

        self.pvtrack = False
        self.non_interacting = False
        self.reverse = True

        self.man = False

    def run_for(self,dt):
        if (self._firstRun):
            self._lastPv = self._nextPv
            self._lastSp = self._nextSp
            self._lastFF = self._nextFF
        dPv = self._nextPv - self._lastPv

        # change in erorr is change in PV minus change in SP
        dErr = dPv - (self._nextSp - self._lastSp)
        dOP = self.ff_gain*(self._nextFF - self._lastFF)/self.opRange

        self._lastPv = self._nextPv
        self._lastSp = self._nextSp
        self._lastFF = self._nextFF

        
        if self.reverse:
            act = -1
        else:
            act = 1

        if not self.man:
            # For derivative action we must calculate the rate of change of the
            # derivative term (pv or error) and then calculate the change in the
            # rate of change
            if (self.D_onErr):
                dXdt = dErr/dt
            else:
                dXdt = dPv/dt
            dDD = dXdt - self._lastDxDt # TODO: Divide by dt?
            self._lastDxDt = dXdt

            # dOP is scaled delta OP. Reverse acting is the norm, that is why
            # everything is negative

            # Integral Action
            if self.Ti <= 0:
                Ki = 0
            elif self.non_interacting:
                Ki = 1/self.Ti
            else:
                Ki = self._K/self.Ti

            dOP += act*dt*Ki*(self._nextPv - self._nextSp)/self.pvRange

            # Proportional Action
            if self.K_onErr:
                dOP += act*self._K*(dErr/self.pvRange)
            else:
                dOP += act*self._K*(dPv/self.pvRange)

            # Derivative Action
            if self.non_interacting:
                Kd = self.Td
            else:
                Kd = self._K* self.Td
            dOP += act*Kd*dDD/self.pvRange

            self.op += dOP*self.opRange
        else:
            if self.pvtrack:
                self.set_sp(self._nextPv,init=True)

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
    tNext        : time for next move

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
        self.tNext = 1.0
        self.opRange = (0.0, 1.0)

    def run_for(self,dt):
        dtHr = dt/3600

        self.tMove += dtHr

        if self.tMove < self.reactionTime[0]:
            return

        if   ( self.pv < self.controlRange[0] ):
            self.op += self.K*(self.controlRange[0] - self.pv)
            self.tMove = 0.0
            self.tNext = numpy.random.uniform()*self.reactionTime[1]
        elif ( self.pv > self.controlRange[1] ):
            self.op += self.K*(self.controlRange[1] - self.pv)
            self.tMove = 0.0
            self.tNext = numpy.random.uniform()*self.reactionTime[1]
        elif ( self.tMove > self.tNext ):
            sp = (self.controlRange[0] + self.controlRange[1])/2
            self.op += 0.5*self.K*(sp - self.pv)
            self.tMove = 0.0
            self.tNext = numpy.random.uniform()*self.reactionTime[1]

        self.op = max(self.opRange[0],min(self.opRange[1],self.op))
        

class ValveCharProcess(AbstractProcess):
    '''
    Simulates effect of a non-linear field actuator.
    Equal percentage valves are installed to make processes that have a
    non-linear characteristic, linear.  Sometimes equal percentage valves are
    installed in the wrong place and the valve introduces non-linearities to the
    system.

    Example: Centrifugal pump with flow controller at outlet.  This is an
    example where you should use an equal percentage valve.  To simulate the
    effect of an installed linear valve use 'quick-opening' characteristic'

    Simulated Parameters:
    ---------------------
    form : string ('linear','equal-percentage','quick-opening')
           lamda function
        Form of non-linearity.  Let l be controller output and f be effective
        output in field.
            linear : f = l
            equal-percentage : f = R^(l-1)
            quick-opening: f = sqrt(l)
            quick-opening-s: f = (l)*(1/s)
        You can also pass a lambda function as the form.  The function should
        convert an input value of range 0 to 1 to an output in range 0 to 1, the
        input will be scaled with opRange before passed to the lambda function.
    R : float
        Characteristing for equal-percentage, usually between 20 and 50 (default 50)
    s : float
        Tunable characteristic for quick-opening-s, if s=2 then it is equivalent
        to quick-opening
    opRange : float
        Output range (default 100)

    Simulated Inputs:
    -----------------
    input : float
        Controller output

    Simulated Outputs:
    ------------------
    out : float
        Valve position
    '''


    def __init__(self,form='linear'):
        super().__init__()
        self.R=50.0
        self.s=2.0
        self.opRange=100
        self.form = form

        self.input = 0.0
        self.out = 0.0

    def run_for(self,dt):
        if callable(self.form):
            self.out = self.opRange * self.form(self.input/self.opRange)
        elif self.form == 'linear':
            self.out = self.input
        elif self.form == 'equal-percentage':
            self.out = self.opRange * self.R**((self.input/self.opRange)-1)
        elif self.form == 'quick-opening':
            self.out = self.opRange * numpy.sqrt(self.input/self.opRange)
        elif self.form == 'quick-opening-s':
            self.out = self.opRange * (self.input/self.opRange)**(1/self.s)
        else:
            raise Exception('Unidentified characteristic: ' + self.form)


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

    @property
    def output(self):
        return self._output
    @output.setter
    def output(self,value):
        self._output = value
        self._input = value
    @property
    def input(self):
        return self._input
    @input.setter
    def input(self,value):
        self._input = value
        self._output = value


    def __init__(self):
        super().__init__()

        self.sigma = 0.05
        self._input = 0.0
        self._output = 0.0

    def run_for(self,dt):
        self._output = numpy.random.normal(self._input,self.sigma)
        pass

class BrownNoiseProcess(AbstractProcess):
    '''
    Simulates process noise. Process noise is Brown noise (I think) (I think
    process noise is brown noise and I think what I did here is brown noise, but
    I won't bet money on it.)

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


class SinNoiseProcess(AbstractProcess):
    '''
    Add a sinusoid noise signal to input
    
    Simulation Parameters:
    ----------------------
    amplitude : Amplitude of sine wave (EU)
    period    : Period if sine wave (seconds)

    Simulation Inputs:
    input     : Process Input (EU)

    Simulation Outputs:
    -------------------
    output    : Output (noise + input)

    '''

    @property
    def output(self):
        return self.input + self.amplitude*numpy.sin(self._x)
    @output.setter
    def output(self,value):
        self.input = value - self.amplitude*numpy.sin(self._x)


    def __init__(self,period=60,amplitude=1):
        super().__init__()
        self.period = period
        self.amplitude = amplitude
        self._x = 0
        self.input = 0

    def run_for(self,dt):
        self._x += (dt/self.period)*2*numpy.pi

        
        
        





class MathAddProcess(AbstractProcess):
    '''
    Adds inputs
    
    Simulation Parameters:
    ----------------------
    input : numpy array
        input values
    scale : numpy array
        scale values for each input

    Simulation Outputs:
    -------------------
    output : sum of inputs
    '''

    @property
    def output(self):
        return sum(self.input * self.scale)
    @output.setter
    def output(self,value):
        for i in range(len(self.input)):
            self.input[i] = value/(len(self.input)*self.scale[i])

    def __init__(self,n_inputs=3):
        '''
        Parameters:
        -----------
        n_inputs: number of inputs
        '''
        super().__init__()

        self.input = numpy.zeros(n_inputs)
        self.scale = numpy.ones(n_inputs)

    def run_for(self,dt):
        pass




class SelectProcess(AbstractProcess):
    '''
    Selects an input as the output based on selection criteria.

    Simulation Parameters:
    ----------------------
    input : list
        Input values
    seltype : string
        Selection Type
        hi     : select minimum of inputs
        lo     : select maximum of inputs
        median : select median, if even number of inputs, select avg of middle
                 two inputs
    rate : Rate at which output is ramped to selected input [EU/sec].

    Simulation States:
    ------------------
    output : Output

    '''

    @property
    def output(self):
        return self._output

    @output.setter
    def output(self,value):
        self._output = value
        for i in range(len(self.input)):
            self.input[i] = value

    def __init__(self,n_inputs=2):
        '''
        Parameters:
        -----------
        n_inputs : numper of inputs
        '''
        super().__init__()
        self.input = [0.0]*n_inputs
        self.seltype = 'hi'
        self.rate = float('inf')

        self._output = 0

    def run_for(self,dt):
        out = None
        if self.seltype == 'hi':
            out = max(self.input)
        elif self.seltype == 'lo':
            out = min(self.input)
        elif self.seltype == 'median':
            out = numpy.median(self.input)
        else:
            raise Exception("Unknown seltype {}".format(self.seltype))


        self._output = min(
                         max(out,  self._output - self.rate*dt),
                         self._output + self.rate*dt
                       )


            
        
class PLTProcess(AbstractProcess):
    '''
    Piecewise Linear Transform
    Input is transformed across piecewise linear function

    Simulation Parameters:
    ----------------------
    ipLimits : (float,float)
               Input limits.  Input is clamped to these limits
    opLimits : (float,float)
               Output limits. Output is clamped to these limits

    Note: See set_transform function for definition of transform

    Simulation Inputs:
    ------------------
    input

    Simulation Outputs:
    -------------------
    output
    '''

    @property
    def output(self):
        x = max(self.ipLimits[0],min(self.ipLimits[1],self.input))
        y = self._transform(x)
        y = max(self.opLimits[0],min(self.opLimits[1],y))
        return self._transform(self.input)

    @output.setter
    def output(self,value):
        y = max(self.opLimits[0],min(self.opLimits[1],value))
        x = self._revTransform(y)
        self.input = max(self.ipLimits[0],min(self.ipLimits[1],x))

    def set_transform(self,X,Y):
        '''
        Set the transform. A scipy.interp.interp1d function is created from
        input- and output transfroms.

        Parameters:
        -----------
        X : list.  Input coordinates. Must be monotonically increasing.
        Y : list.  Output coordinates. Initialisation of output is unpredictable
                   if not monotonically increasing order.
        '''

        self._transform = scipy.interpolate.interp1d(
            x=X,
            y=Y,
            fill_value='extrapolate',
            assume_sorted=True,
        )
        self._revTransform = scipy.interpolate.interp1d(
            x=Y,
            y=X,
            fill_value='extrapolate',
            assume_sorted=True
        )

        
    def __init__(self):
        super().__init__()

        self.input = 50

        self.ipLimits = (-float('inf'),float('inf'))
        self.opLimits = (-float('inf'),float('inf'))

        self.set_transform([0,100],[0,100])


    def run_for(self,dt):
        pass




