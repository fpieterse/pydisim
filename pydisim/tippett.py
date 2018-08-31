'''
This file contains a CSTR and Separator process used in
Distributed Model Predictive Control Based on Dissipativity
Tippett, Bao, 2012

The following discrepancies are noted:
1) The dynamics of the process appears to be much faster than in (Tippett,2012)
2) Tippett (2012) does not give values for all parameters, specifically
flowrates.
3) The composition of the top product in the separator is a function of relative
volatilities. We expect the heat input to affect the top product flowrate but
none of the equations mention this relationship.
4) Tippett (2012) does not define Fp, if it is the bottoms flowrate then we
don't agree with equations 77 and 78.
5) We decided to make the separator overhead flowrate a manipulated variable (F5
in Tippett (2012)). The heat input to the separator affects the temperature but
not much else. It could probably be ignored.
6) Equation 79 in Tippett (2012) needs a relationship to the heat of evaporation
of the material.

'''

import math
from .processes import AbstractProcess, MainProcess, PIDProcess

class CSTRProcess(AbstractProcess):
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
    Fout_F : float
        Output flowrate
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
    @property
    def Fout_F(self):
        return self.F1_F + self.F2_F

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

        self.k1 = 9.97e6
        self.E1 = 5e4
        self.H1 = -6e4
        self.k2 = 9.0e6
        self.E2 = 6e4
        self.H2 = -7e4

        self.Cp = 4200 
        self.J_Cp = 4200 
        self.rho = 1000
        self.J_rho = 1000

        self._xA = 1.0
        self._xB = 0.0

    def run_for(self,dt):
        R = 8.314 # J/molK
        #extent of reactions:

        rX1 = self.k1 * math.exp(-self.E1/(R*(self.T+273.15))) * self.xA
        rX2 = self.k2 * math.exp(-self.E2/(R*(self.T+273.15))) * self.xB

        # rate of change of component A
        dxA = ((self.F1_F/self.V)*(self.F1_xA-self.xA) +
               (self.F2_F/self.V)*(self.F2_xA-self.xA) -
                rX1 )

        dxB = ((self.F1_F/self.V)*(self.F1_xB-self.xB) +
               (self.F2_F/self.V)*(self.F2_xB-self.xB) +
                rX1 - rX2 )

        dT = ( (self.F1_F/self.V)*(self.F1_T-self.T) +
               (self.F2_F/self.V)*(self.F2_T-self.T) -
               (self.H1/self.Cp)*rX1 -
               (self.H2/self.Cp)*rX2 +
                self.UA*3600 * (self.J_T-self.T) / (self.rho*self.Cp*self.V) )

        dTJ = ( ( self.FJ_F / self.J_V )*( self.FJ_T - self.J_T ) -
                self.UA*3600 * (self.J_T-self.T) / (self.J_rho*self.J_Cp*self.J_V) )

       
        dthr = dt/3600
        self._xA += dxA*dthr
        self._xB += dxB*dthr
        self.T += dT*dthr
        self.J_T += dTJ*dthr

class SeparatorProcess(AbstractProcess):
    '''
    Separator 

    Simualtion Paramaters:
    -----------
    V : float
        Volume of separator
    rv_[A|B|C] : float
        Relative Volatility of species
    Cp : float
        Heat capacity (J/kg/K)
    rho : float
        Density (kg/m3)


    Simulation Inputs:
    ------------------
    Q : float
        Heat input to separator (W)
    Fin_F : float
        Feed flowrate
    Fin_T : float
        Feed temperature
    Fin_x[A|B] : float
        Feed composition
    Ftop_F : float
        Top product flowrate


    Simulated States:
    -----------------
    x[A|B] : float
        Fractions of components A, B, C. Can only be set using the
        set_components function.
    T : float
        Temperature (°C)

    Simulated Outputs:
    ------------------
    xC : float
        Fraction of component C
    Ftop_x[A|B|C] : float
        Top product fractions
    Fbot_F : float
        Bottom product flowrate
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
    @property
    def Ftop_xC(self):
        return 1 - self.Ftop_xA - self.Ftop_xB
    @property
    def Fbot_F(self):
        return self.Fin_F - self.Ftop_F

    def set_components(self,xA,xB):
        if (xA + xB) > 1.0:
            raise ValueError("Fractions of A and B bigger than 1.0")
        self._xA = xA
        self._xB = xB

    def __init__(self):
        super().__init__()

        self.V = 1.0
        self.rv_A = 3.5
        self.rv_B = 1.0
        self.rv_C = 0.5
        self.Cp = 4200
        self.rho = 1000

        self._xA = 0.4
        self._xB = 0.4
        self.T = 25

        self.Fin_F = 0.0
        self.Fin_T = 25
        self.Fin_xA = 0.4
        self.Fin_xB = 0.4
        self.Ftop_F = 0.0
        self.Q = 0.0

        self._calculate_top_fractions()

    def _calculate_top_fractions(self):
        rxa = self.rv_A*self._xA
        rxb = self.rv_B*self._xB
        rxc = self.rv_C*self.xC
        rsum = rxa+rxb+rxc
        self.Ftop_xA = rxa/rsum
        self.Ftop_xB = rxb/rsum

    def run_for(self,dt):
        self._calculate_top_fractions()

        dxA = ( (self.Fin_F/self.V)*(self.Fin_xA - self._xA) -
                (self.Ftop_F/self.V)*(self.Ftop_xA - self._xA) )

        dxB = ( (self.Fin_F/self.V)*(self.Fin_xB - self._xB) -
                (self.Ftop_F/self.V)*(self.Ftop_xB - self._xB) )

        dT = ( (self.Fin_F/self.V)*(self.Fin_T - self.T) +
               (self.Q*3600)/(self.rho*self.Cp*self.V) )

        dthr = dt/3600

        self._xA += dxA*dthr
        self._xB += dxB*dthr
        self.T += dT*dthr

def CreateProcess():
    '''
    Creates the process described in Tippett (2012).

    Returns a MainProcess instance with 2 CSTRs in series and a separator. PID
    controllers are included to ensure stability.
    '''

    rx1 = CSTRProcess()
    rx1.V = 1.0
    rx1.J_V = 0.1
    rx1.UA = 2509.8 

    rx1.F1_F = 5.0
    rx1.F1_T = 150.0
    rx1.FJ_F = 0.5
    rx1.FJ_T = 280
    rx1.J_T = 190
    rx1.T = 170

    rx1.set_components(0.38,0.56)


    rx2 = CSTRProcess()
    rx2.V = 0.5
    rx2.J_V = 0.05
    rx2.UA = 2774

    rx2.F1_F = 1.0
    rx2.F1_T = 150.0
    rx2.FJ_F = 0.2
    rx2.FJ_T = 280
    rx2.J_T = 175
    rx2.T = 170

    rx2.set_components(0.31,0.63)


    sep = SeparatorProcess()
    sep.V = 1.0
    sep.Q = 20000
    sep.T = 178

    sep.Ftop_F = 10.0
    rx1.F2_F = 10.0
    rx1.F2_T = 178

    sep.set_components(0.17,0.75)



    # Separator bottoms B-Spec controller

    xc02 = PIDProcess()
    xc02.sp = 0.75
    xc02.pvRange = 0.2
    xc02.opRange = 10.0
    xc02.opLimits = (0.0,10.0)
    xc02.op = rx1.F1_F
    xc02.K = -1.0
    xc02.Ti = 800

    tc01 = PIDProcess()
    tc01.sp = 170
    tc01.op = 0.5
    tc01.opLimits = (0,2.0)
    tc01.opRange = 2.0
    tc01.pvRange = 20.0
    tc01.K = 2.0
    tc01.Ti = 1500

    tc02 = PIDProcess()
    tc02.sp = 170
    tc02.op = 0.2
    tc02.opLimits = (0,1.0)
    tc02.opRange = 1.0
    tc02.pvRange = 20.0


    main = MainProcess()

    main.AddProcess(rx1)
    main.AddProcess(tc01)
    main.AddProcess(rx2)
    main.AddProcess(tc02)

    main.AddProcess(sep)
    main.AddProcess(xc02)


    main.AddConnection(rx1,'Fout_F',rx2,'F2_F')
    main.AddConnection(rx1,'T',rx2,'F2_T')
    main.AddConnection(rx1,'xA',rx2,'F2_xA')
    main.AddConnection(rx1,'xB',rx2,'F2_xB')

    main.AddConnection(rx1,'T',tc01,'pv')
    main.AddConnection(tc01,'op',rx1,'FJ_F')

    main.AddConnection(rx2,'Fout_F',sep,'Fin_F')
    main.AddConnection(rx2,'T',sep,'Fin_T')
    main.AddConnection(rx2,'xA',sep,'Fin_xA')
    main.AddConnection(rx2,'xB',sep,'Fin_xB')

    main.AddConnection(rx2,'T',tc02,'pv')
    main.AddConnection(tc02,'op',rx2,'FJ_F')
    
    main.AddConnection(sep,'Ftop_F',rx1,'F2_F')
    main.AddConnection(sep,'T',rx1,'F2_T')
    main.AddConnection(sep,'Ftop_xA',rx1,'F2_xA')
    main.AddConnection(sep,'Ftop_xB',rx1,'F2_xB')

    main.AddConnection(sep,'xB',xc02,'pv')
    main.AddConnection(xc02,'op',rx1,'F1_F')

    main.rx1 = rx1
    main.rx2 = rx2
    main.sep = sep
    main.xc02 = xc02
    main.tc01 = tc01
    main.tc02 = tc02

    for i in range(1000):
        main.run_for(10)

    return main

