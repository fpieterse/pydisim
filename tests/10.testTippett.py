#!/usr/bin/env python
'''
Tests for tippett processes
'''

import sys
sys.path.append('../')

import matplotlib.pyplot as plt

from pydisim.tippett import *
from pydisim.tools import Recorder

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

# Tippett Process >>>

mp = CreateProcess()

t = 0
dt = 10

time = []
rx1CompRecorder = Recorder(['xA_1','xB_1','xC_1'])
rx2CompRecorder = Recorder(['xA_2','xB_2','xC_2'])
tempRecorder = Recorder(['T1','T2','T3','TJ1','TJ2'])


for i in range(100):
    t += dt
    mp.RunFor(dt)
    time.append(t)
    rx1CompRecorder.record(mp.rx1.xA, mp.rx1.xB, mp.rx1.xC)
    rx2CompRecorder.record(mp.rx2.xA, mp.rx2.xB, mp.rx2.xC)

    tempRecorder.record(mp.rx1.T, mp.rx2.T, mp.sep.T, mp.rx1.J_T, mp.rx2.J_T)


fig,axes = plt.subplots(2,1,sharex=True)
rx1CompRecorder.plot(axes[0],time)
rx2CompRecorder.plot(axes[1],time)
fig.tight_layout()

fig2 = plt.figure(2)
tempRecorder.plot(plt,time)

plt.show()

# <<<

# doSep >>>
if doSep:
    sep = SeparatorProcess()
    sep.T = 230
    sep.Fin_F = 1.0
    sep.Fin_T = 230
    sep.Ftop_F = 0.5
    sep.Q = 2000

    time = []
    T = []
    xA = []
    xB = []
    xC = []
    xAtop = []
    xBtop = []
    xCtop = []

    t = 0
    dt = 5
    for i in range(1000):
        t += dt
        time.append(t)
        sep.RunFor(dt)
        T.append(sep.T)
        xA.append(sep.xA)
        xB.append(sep.xB)
        xC.append(sep.xC)
        xAtop.append(sep.Ftop_xA)
        xBtop.append(sep.Ftop_xB)
        xCtop.append(sep.Ftop_xC)

    fig, axes = plt.subplots(2,1,sharex=True)
    axes[0].plot(time,T)
    axes[1].plot(time,xA,'b-',label='xA')
    axes[1].plot(time,xB,'g-',label='xB')
    axes[1].plot(time,xC,'r-',label='xC')
    axes[1].plot(time,xAtop,'b--',label='xAtop')
    axes[1].plot(time,xBtop,'g--',label='xBtop')
    axes[1].plot(time,xCtop,'r--',label='xCtop')
    axes[1].legend()

    fig.tight_layout()

    plt.show()
#<<<

# doCSTR >>>
if doCSTR:
    #---- Reactors -----------------------------------------------------------
    # CSTR with no inflow
    rx = CSTRProcess()
    rx.T = 230
    rx.J_T = 300
    rx.FJ_T = 300
    rx.FJ_F = 1.0

    time = []
    T = []
    TJ = []
    xA = []
    xB = []
    xC = []

    dt = 5
    t = 0
    for i in range(3600):
        t += dt
        time.append(t)
        rx.RunFor(dt)
        T.append(rx.T)
        TJ.append(rx.J_T)
        xA.append(rx.xA)
        xB.append(rx.xB)
        xC.append(rx.xC)

    plt.plot(time,xA,label='xA')
    plt.plot(time,xB,label='xB')
    plt.plot(time,xC,label='xC')
    plt.legend()
    plt.show()

    plt.plot(time,T,label='T')
    plt.plot(time,TJ,label='Tjacket')
    plt.legend()
    plt.show()

# <<<

# vim:fdm=marker:fmr=>>>,<<<:
