
class AbstractProcess:
    '''
    Abstract process unit.

    Represents a process unit with inputs and outputs that can be simulated for
    a given number of seconds by the MainProcess process unit (that is itself
    also a process unit).

    Inherited classes must implement the RunFor-function.
    Inherited classes must execute the super-class' __init__-function that
    creates a list of connection strings used by the MainProcess.

    '''

    def __init__(self):
        # Record of connection info used by MainProcess
        # This will be a list of tuples of the form
        # (outputName, downstreamProcess, downstreamInputName)
        self._connectionInfo = []

    def RunFor(self,dt):
        '''
        Run process for dt seconds. Must be overloaded
        '''
        raise NotImplementedError

class MainProcess(AbstractProcess):
    '''
    The main process represents a collection of subprocesses and their
    connections and relationships
    '''
    def __init__(self):
        super().__init__()
        self._subProcesses = []

    def AddProcess(self,process):
        '''
        Adds a process to the list of processes
        '''
        self._subProcesses.append(process)

    def AddConnection(self,upstreamProcess,
                           upstreamProcessOutput,
                           downstreamProcess,
                           downstreamProcessInput):
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
        '''
        upstreamProcess._connectionInfo.append(
            (upstreamProcessOutput,
             downstreamProcess,
             downstreamProcessInput))

    def RunFor(self,dt):
        '''
        Run main process for dt seconds.
        This involves running each subprocess (in the order they were added to
        the main process) and then writing that subprocess' outputs to the
        relevant downstream processes.
        '''
        for proc in self._subProcesses:
            proc.RunFor(dt)
            for c in proc._connectionInfo:
                setattr(c[1],c[2],getattr(proc,c[0]))


class HoldupProcess(AbstractProcess):
    '''
    Simulated tank.
    Tank level can go beyon 0-100%


    Creation:
    ---------
    TankProcess(vol,level_0)
    vol : float
        volume of vessel (m3)
    level_0 : float
        starting level (%)

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

    Simulated Outputs:
    ------------------
    level : float
        liquid level in tank (%)
    '''

    def __init__(self,vol,level_0):
        super().__init__()

        self._vol = vol
        self.level = level_0
        self.cVol = self._vol*self.level/100.0

        self.fIn = 0
        self.fOut = 0

    def RunFor(self,dt):
        self.cVol += dt*(self.fIn - self.fOut)/3600.0
        self.level = 100.0 * self.cVol/self._vol

import math
class CSTRProcessI(AbstractProcess):
    '''
    Constant holdup CSTR process.

    Reactor has two feeds that are mixtures of A, B, and C. The reaction in the
    reactor is A -> B -> C and B is the valuable product. The reactor
    temperature is controlled using a jacket. The flowrate of material through
    the jacket can be adjusted.

    Parameters:
    -----------
    V : float
        Volume of reactor
    J_V : float
        Volume of jacket
    k1 : float
        Reaction rate for A -> B (1/hr)
    E1 : float
        Activation eneregy for A -> B (J/mol)
    H1 : float
        Heat of reaction (J/mol)
    k2 : float
        Reaction rate for B -> C (1/hr)
    E2 : float
        Activation eneregy for B -> C (J/mol)
    H2 : float
        Activation energy (J/mol)
    UA : float
        Jacket heat transfer coefficient (W/K)
    Cp : float
        Heat capacity in the reactor (J/kg/K)
    J_Cp : float
        Heat capacity of jacket fluid (J/kg/K)
    rho : float
        Material density
    J_rho : float
        Jacket material density


    Simulated Inputs:
    -----------------
    Fn_F : float
        Fn flow rate (m3/h)
    Fn_T : float
        Fn flow temperature (degC)
    Fn_x[A|B] : float
        Fractions of components A and B in stream n


    Simulated States:
    -----------------
    T : float
        Temperature (degC)
    J_T : float
        Jacket temperature (degC)
    x[A|B] : float
        Fraction of components A and B.
        Fractions can only be set using the set_components function


    Simpulated Outputs:
    -------------------
    xC : float
        Fraction of component C
    '''

    @property
    def xA(self):
        return self._xA
    @property
    def xB(self):
        return self._xB
    @property
    def xC(self):
        return 1 - self._xA - self._xB

    def set_components(self,xA,xB):
        if (xA + xB) > 1.0:
            raise ValueError("Fractions of A and B bigger than 1.0")
        self._xA = xA
        self._xB = xB

    def __init__(self):
        super().__init__()

        self.V = 1.0
        self.J_V = 0.1
        self.UA = 2509.8

        self.T = 25.0
        self.J_T = 25.0

        self.F1_F = 0.0
        self.F1_T = 25.0
        self.F1_xA = 1.0
        self.F1_xB = 0.0

        self.F2_F = 0.0
        self.F2_T = 25.0
        self.F2_xA = 1.0
        self.F2_xB = 0.0

        self.FJ_F = 0.0
        self.FJ_T = 25.0

        self.k1 = 9.97e5
        self.E1 = 5e4
        self.H1 = -6e4
        self.k2 = 9.0e5
        self.E2 = 6e4
        self.H2 = -7e4

        self.Cp = 4200 
        self.J_Cp = 4200 
        self.rho = 1000
        self.J_rho = 1000

        self._xA = 1.0
        self._xB = 0.0

    def RunFor(self,dt):
        R = 8.314 # J/molK
        #extent of reactions:

        print(-self.E1/(R*(self.T+273.15)))

        rX1 = self.k1 * math.exp(-self.E1/(R*(self.T+273.15))) * self.xA
        rX2 = self.k2 * math.exp(-self.E2/(R*(self.T+273.15))) * self.xB
        print("rX1",rX1,"rX2",rX2)

        # rate of change of component A
        dxA = ((self.F1_F/self.V)*(self.F1_xA-self.xA) +
               (self.F2_F/self.V)*(self.F2_xA-self.xA) -
                rX1 )

        dxB = ((self.F1_F/self.V)*(self.F1_xB-self.xB) +
               (self.F2_F/self.V)*(self.F2_xB-self.xB) +
                rX1 - rX2 )

        dT = ( ( self.F1_F / self.V )*( self.F1_T - self.T ) +
               ( self.F2_F / self.V )*( self.F2_T - self.T ) -
               ( self.H1 / self.Cp )*rX1 -
               ( self.H2 / self.Cp )*rX2 +
               self.UA * ( self.J_T - self.T )/( self.rho * self.Cp * self.V ) )

        dTJ = ( ( self.FJ_F / self.J_V )*( self.FJ_T - self.J_T )/
                ( self.J_rho * self.J_Cp * self.J_V) )

       
        dthr = dt/3600
        print("xA",self.xA,self._xA)
        self._xA += dxA*dthr
        print("xA",self.xA,self._xA)
        self._xB += dxB*dthr
        self.T += dT*dthr
        self.J_T += dTJ*dthr


class PIDProcess(AbstractProcess):
    '''
    PID controller

    Creation:
    ---------
    PIDProcess()

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


    def RunFor(self,dt):
        if (self._firstRun):
            self._lastPv = self._nextPv
            self._lastSp = self._nextSp
        dPv = self._nextPv - self._lastPv

        # change in erorr is change in PV minus change in SP
        dErr = dPv - (self._nextSp - self._lastSp)

        self._lastPv = self._nextPv
        self._lastSp = self._nextSp
       
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
        

