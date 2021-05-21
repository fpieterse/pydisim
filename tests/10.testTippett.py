#!/usr/bin/env python
'''
Tests for tippett processes
'''

import sys
sys.path.append('../')

import matplotlib.pyplot as plt

from pydisim.tippett import *
from pydisim.tools import *

doCSTR = False
doSep = False
for a in sys.argv[1:]:
    if a == 'CSTR':
        print("Testing CSTR")
        doCSTR = True
    elif a == 'Sep':
        print("Testing Sep")
        doSep = True
    else:
        print("Unknown argument " + a)

# doCSTR >>>
if doCSTR:
    #---- Reactors -----------------------------------------------------------
    # CSTR with no inflow
    pm = ProcessManager()
    rx = CSTRProcess()
    rx.FJ_F = 1.0
    rx.FJ_T = 300
    rx.F1_F = 1.0
    rx.set_components(0.14,0.65)
    rx.T = 146
    rx.J_T = 194

    rec_CSTR = RecorderProcess(['T','TJ','xA','xB','xC','F1'],rec_int=10)
    rx.add_connection('T',rec_CSTR,'T')
    rx.add_connection('J_T',rec_CSTR,'TJ')
    rx.add_connection('xA',rec_CSTR,'xA')
    rx.add_connection('xB',rec_CSTR,'xB')
    rx.add_connection('xC',rec_CSTR,'xC')
    rx.add_connection('F1_F',rec_CSTR,'F1')
    pm.exec_int = 1

    pm.run_process(30)
    rx.F1_F = 1.2
    pm.run_process(30)
    rx.FJ_F = 1.2
    pm.run_process(60)

    rec_CSTR.plot(subplots=True)

    print('''
The reaction should start close to a steady state.
After that the feed will increase which will cause the concentration of A to
increase and the temperature to drop.  The jacket flowrate then increases which
will increase the temperature and the reaction rate again.
''')
    plt.show()

    exit()
# <<<

# doSep >>>
if doSep:
    pm = ProcessManager()
    sep = SeparatorProcess()
    sep.T = 70
    sep.Fin_F = 1.0
    sep.Fin_T = 90
    sep.Ftop_F = 0.5
    sep.Q = 200

    rec_Sep = RecorderProcess(['T','Ftop','Fbot','xA','xB','Q'],rec_int=10)
    rec_SepX = RecorderProcess(['xA',
                                'xB',
                                'xC',
                                'xA_top',
                                'xB_top',
                                'xC_top'],rec_int=10)

    sep.add_connection('T',rec_Sep,'T')
    sep.add_connection('Ftop_F',rec_Sep,'Ftop')
    sep.add_connection('Fbot_F',rec_Sep,'Fbot')
    sep.add_connection('xA',rec_Sep,'xA')
    sep.add_connection('xB',rec_Sep,'xB')
    sep.add_connection('Q',rec_Sep,'Q')

    sep.add_connection('xA',rec_SepX,'xA')
    sep.add_connection('xB',rec_SepX,'xB')
    sep.add_connection('xC',rec_SepX,'xC')
    sep.add_connection('Ftop_xA',rec_SepX,'xA_top')
    sep.add_connection('Ftop_xB',rec_SepX,'xB_top')
    sep.add_connection('Ftop_xC',rec_SepX,'xC_top')

    pm.exec_int = 10

    pm.run_process(60)
    sep.Q = 500
    pm.run_process(20)
    sep.Q = 200
    pm.run_process(120)

    rec_Sep.plot(subplots=True)

    rec_SepX.plot()


    print('''The temperature should start below the boiling point of 100Â°C and
increase to the boiling point.  At this point the top flow will increase and
the composition in the separator will change to lower xA and higher xC.''')
    plt.show()

    exit()
#<<<

# Tippett Process >>>

pm, tpt_proc = CreateProcess()

rx1CompRecorder = RecorderProcess(['xA_1','xB_1','xC_1'],rec_int=10)
tpt_proc.rx1.add_connection('xA',rx1CompRecorder,'xA_1')
tpt_proc.rx1.add_connection('xB',rx1CompRecorder,'xB_1')
tpt_proc.rx1.add_connection('xC',rx1CompRecorder,'xC_1')

rx2CompRecorder = RecorderProcess(['xA_2','xB_2','xC_2'],rec_int=10)
tpt_proc.rx2.add_connection('xA',rx2CompRecorder,'xA_2')
tpt_proc.rx2.add_connection('xB',rx2CompRecorder,'xB_2')
tpt_proc.rx2.add_connection('xC',rx2CompRecorder,'xC_2')

sepTCompRecorder = RecorderProcess(['Top xA','Top xB','Top xC'],rec_int=10)
tpt_proc.sep.add_connection('Ftop_xA',sepTCompRecorder,'Top xA')
tpt_proc.sep.add_connection('Ftop_xB',sepTCompRecorder,'Top xB')
tpt_proc.sep.add_connection('Ftop_xC',sepTCompRecorder,'Top xC')

sepBCompRecorder = RecorderProcess(['Bot xA','Bot xB','Bot xC'],rec_int=10)
tpt_proc.sep.add_connection('xA',sepBCompRecorder,'Bot xA')
tpt_proc.sep.add_connection('xB',sepBCompRecorder,'Bot xB')
tpt_proc.sep.add_connection('xC',sepBCompRecorder,'Bot xC')

tempRecorder = RecorderProcess(['T1','T2','T3','TJ1','TJ2'],rec_int=10)
tpt_proc.rx1.add_connection('T',tempRecorder,'T1')
tpt_proc.rx2.add_connection('T',tempRecorder,'T2')
tpt_proc.sep.add_connection('T',tempRecorder,'T3')
tpt_proc.rx1.add_connection('J_T',tempRecorder,'TJ1')
tpt_proc.rx2.add_connection('J_T',tempRecorder,'TJ2')

tc01Rec = PIDRecorderProcess('tc01',tpt_proc.tc01,10)
tc02Rec = PIDRecorderProcess('tc02',tpt_proc.tc02,10)

pm.run_process(60)

tpt_proc.tc01.sp += 4
pm.run_process(60)
tpt_proc.tc02.sp += 4
pm.run_process(60)

fig,ax = plt.subplots(nrows=2,sharex=True)
tc01Rec.plot(ax[0])
tc02Rec.plot(ax[1])



fig,axes = plt.subplots(2,2,sharex=True)
rx1CompRecorder.plot(axes[0,0])
rx2CompRecorder.plot(axes[1,0])
sepTCompRecorder.plot(axes[0,1])
sepBCompRecorder.plot(axes[1,1])
fig.tight_layout()

fig2 = plt.figure(2)
tempRecorder.plot()


plt.show()

# <<<



# vim:fdm=marker:fmr=>>>,<<<:fileencoding=utf-8:
