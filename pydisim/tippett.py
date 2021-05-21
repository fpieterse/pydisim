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
from .processes import AbstractProcess, \
                       PIDProcess, \
                       ProcessManager, \
                       GaussNoiseProcess

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
        Heat of reaction (kJ/mol)
    k2 : float
        Reaction rate for B -> C (1/hr)
    E2 : float
        Activation eneregy for B -> C (J/mol)
    H2 : float
        Activation energy (kJ/mol)
    UA : float
        Jacket heat transfer coefficient (kW/K)
    Cp : float
        Heat capacity in the reactor (kJ/kg/K)
    J_Cp : float
        Heat capacity of jacket fluid (kJ/kg/K)
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
        self.UA = 2.5098

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
        self.H1 = -6e1
        self.k2 = 9.0e6
        self.E2 = 6e4
        self.H2 = -7e1

        self.Cp = 4.200
        self.J_Cp = 4.200
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
    rho : float
        Density (kg/m3)
    Cp : float
        Heat Capacity (kJ/kgK)
    Tbp : float
        Boiling Point (ÃÂ°C)
    Hvap : float
        Heat of Vaporisation (kJ/kg)


    Simulation Inputs:
    ------------------
    Q : float
        Heat input to separator (kW)
    Fin_F : float
        Feed flowrate (m3/hr)
    Fin_T : float
        Feed temperature
    Fin_x[A|B] : float
        Feed composition


    Simulated States:
    -----------------
    x[A|B] : float
        Fractions of components A, B, C. Can only be set using the
        set_components function.
    T : float
        Temperature (ÃÂ°C)
        

    Simulated Outputs:
    ------------------
    xC : float
        Fraction of component C
    Ftop_x[A|B|C] : float
        Top product fractions
    Ftop_F : float
        Top product flowrate (m3/h)
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
        self.Cp = 4.2
        self.rho = 1000
        self.Hvap = 2200

        self._xA = 0.4
        self._xB = 0.4
        self.T = 25
        self.Tbp = 100

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

        dthr = dt/3600

        dT = ( (self.Fin_F/self.V)*(self.Fin_T - self.T) +
               (self.Q*3600)/(self.rho*self.Cp*self.V) )

        # how much extra temperature above boiling point do we have
        t_extra = dT*dthr - (self.Tbp - self.T)
        if t_extra > 0:
            self.T = self.Tbp
            q_extra = t_extra*self.rho*self.Cp*self.V
            self.Ftop_F = (q_extra/(self.Hvap*self.rho) ) /dthr
            if self.Ftop_F > self.Fin_F:
                self.Ftop_F = self.Fin_F
        else:
            self.T += dT*dthr
            self.Ftop_F = 0

        dxA = ( (self.Fin_F/self.V)*(self.Fin_xA - self._xA) -
                (self.Ftop_F/self.V)*(self.Ftop_xA - self._xA) )

        dxB = ( (self.Fin_F/self.V)*(self.Fin_xB - self._xB) -
                (self.Ftop_F/self.V)*(self.Ftop_xB - self._xB) )


        self._xA += dxA*dthr
        self._xB += dxB*dthr

class TippettProcess(AbstractProcess):
    '''
    2 CSTR 1 Seperator process described by Tippett 2012.


    Simulated Inputs:
    -----------------
    F1 : float
        Fresh feed to rx1 (m3/h)
    F2 : float
        Fresh feed to rx2 (m3/h)
    Q : float
        Heat input to separator (kW)


    Simulated Outputs:
    ------------------
    F3 : float
        Feed from rx1 to rx2 (m3/h)
    F4 : float
        Feed to separator (m3/h)
    F5 : float
        Recycle flow (m3/h)
    F6 : float
        Product flow (m3/h)
    x[A|B|C]_prod : float
        Composition of product
    Fprod : float
        Flowrate of component B in product (m3/h)



    SubProcesses:
    -----------
    rx1 : tippett.CSTRProcess
        CSTR 1
    rx2 : tippett.CSTRProcess
        CSTR 2
    sep : tippett.SepProcess
        Seperator
    tc01 : PIDProcess
        CSTR1 temperature controller
    tc02 : PIDProcess
        CSTR2 temeperature controller
        
    '''

    @property
    def F1(self):
        return self.rx1.F1_F
    @F1.setter
    def F1(self,value):
        self.rx1.F1_F = value

    @property
    def F2(self):
        return self.rx2.F1_F
    @F2.setter
    def F2(self,value):
        self.rx2.F1_F = value

    @property
    def Q(self):
        return self.sep.Q
    @Q.setter
    def Q(self,value):
        self.sep.Q = value

    @property
    def F3(self):
        return self.rx2.F2_F
    @property
    def F4(self):
        return self.sep.Fin_F
    @property
    def F5(self):
        return self.sep.Ftop_F
    @property
    def F6(self):
        return self.sep.Fbot_F

    @property
    def xA_prod(self):
        return self.sep.xA
    @property
    def xB_prod(self):
        return self.sep.xB
    @property
    def xC_prod(self):
        return self.sep.xC
    @property
    def Fprod(self):
        return self.sep.xB * self.sep.Fbot_F


    def __init__(self):
        super().__init__()


        # RX1

        self.rx1 = CSTRProcess()
        self.rx1.V = 1.0
        self.rx1.J_V = 0.1
        self.rx1.UA = 5.0

        self.rx1.F1_F = 5.0
        self.rx1.F1_T = 150.0
        self.rx1.FJ_F = 0.5
        self.rx1.FJ_T = 250
        self.rx1.J_T = 190
        self.rx1.T = 170

        self.rx1.F2_F = 10.0
        self.rx1.F2_T = 170

        self.rx1.set_components(0.38,0.56)

        # RX2

        self.rx2 = CSTRProcess()
        self.rx2.V = 0.5
        self.rx2.J_V = 0.05
        self.rx2.UA = 3.0

        self.rx2.F1_F = 1.0
        self.rx2.F1_T = 150.0
        self.rx2.FJ_F = 0.5
        self.rx2.FJ_T = 250
        self.rx2.J_T = 175
        self.rx2.T = 170

        self.rx2.set_components(0.31,0.63)

        self.rx1.add_connection('Fout_F',self.rx2,'F2_F')
        self.rx1.add_connection('T',self.rx2,'F2_T')
        self.rx1.add_connection('xA',self.rx2,'F2_xA')
        self.rx1.add_connection('xB',self.rx2,'F2_xB')

        # SEP

        self.sep = SeparatorProcess()
        self.sep.V = 1.0
        self.sep.Q = 1000

        self.sep.set_components(0.17,0.75)

        self.rx2.add_connection('Fout_F',self.sep,'Fin_F')
        self.rx2.add_connection('T',self.sep,'Fin_T')
        self.rx2.add_connection('xA',self.sep,'Fin_xA')
        self.rx2.add_connection('xB',self.sep,'Fin_xB')

        self.sep.add_connection('Ftop_F',self.rx1,'F2_F')
        self.sep.add_connection('T',self.rx1,'F2_T')
        self.sep.add_connection('Ftop_xA',self.rx1,'F2_xA')
        self.sep.add_connection('Ftop_xB',self.rx1,'F2_xB')

        # TC01 : RX1 temperature controller
        self.tc01 = PIDProcess()
        self.tc01.sp = 155
        self.tc01.opLimits = (0,10.0)
        self.tc01.opRange = 10.0
        self.tc01.pvRange = 20.0
        self.tc01.K = 1.0
        self.tc01.Ti = 400

        self.rx1.add_connection('T',self.tc01,'pv')
        self.tc01.add_connection('op',self.rx1,'FJ_F','<')

        # TC01 : RX2 temperature controller

        self.tc02 = PIDProcess()
        self.tc02.sp = 170
        self.tc02.opLimits = (0,10.0)
        self.tc02.opRange = 10.0
        self.tc02.pvRange = 20.0
        self.tc02.K = 2.0
        self.tc02.Ti = 300

        self.rx2.add_connection('T',self.tc02,'pv')
        self.tc02.add_connection('op',self.rx2,'FJ_F','<')


    def run_for(self,dt):
        # Process Manager will run all subprocess.
        # If subprocess inputs are exposed, they will be connected here.
        pass


def CreateProcess():
    '''
    Creates the process described in Tippett (2012).

    Returns the ProcessManager and the TippettProcess
    '''

    pm = ProcessManager(exec_int=5)
    tp = TippettProcess()
    pm.run_process(120)

    return pm, tp

# vim: fileencoding=utf-8
