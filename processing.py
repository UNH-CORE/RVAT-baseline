"""
This script imports and processes turbine performance and wake data from
March 2013 experiments with the VAT

"""
from __future__ import division, print_function
import pyTDMS
import xlrd
import timeseries
from timeseries import *
import numpy as np
import csv
import matplotlib.pyplot as plt
import time
import matplotlib
from scipy.signal import decimate
from scipy.interpolate import interp1d
import fdiff
import sys
import os

# Some constants
R = 0.5
D = 2*R
H = 1.0
U = 1.0
d_shaft = 0.095
A_t = 1.0
rho = 1000.0
nu = 1e-6
    
def import_testplan():
    wb = xlrd.open_workbook("Test plan, 2013.03 VAT, Rev4.xlsx")
    ws = wb.sheet_by_index(0)
    runs = ws.col_values(0)
    # Find row with run 1 in it
    firstrow = runs.index(1)
    runs = [int(run) for run in runs[firstrow:]]
    tsr = ws.col_values(1)[firstrow:]
    y_R = ws.col_values(2)[firstrow:]
    z_H = ws.col_values(4)[firstrow:]
    return {"runs" : np.asarray(runs), "tsr" : np.asarray(tsr), 
            "y/R" : np.asarray(y_R), "z/H" : np.asarray(z_H)}
    
def loadvec(run):
    data = np.loadtxt("Vectrino/vec" + str(run) + ".dat")
    t = data[:,0]
    u = data[:,3]
    v = data[:,4]
    w = data[:,6]
    return t, u, v, w

def loadvectemp(run):
    with open("Vectrino/vec" + str(run) + ".hdr") as f:
        temp = f.readlines()[117].split()[1]
    return float(temp)
                          
def loadtdms(run):
    filename = "TDMS/run" + str(run) + ".tdms"
    objects, rawdata = pyTDMS.read(filename)
    Ttrans = np.asarray(rawdata[b"/'Untitled'/'TorqueTrans'"])
    Tarm = np.asarray(rawdata[b"/'Untitled'/'TorqueArm'"])
    dragL = np.asarray(rawdata[b"/'Untitled'/'DragL'"])
    dragR = np.asarray(rawdata[b"/'Untitled'/'DragR'"])
    Tarm = Tarm - np.mean(Tarm[0:2000])
    tareDrag = 49.2407
    drag = dragL + dragR
    drag = drag - np.mean(drag[0:2000])
    drag = drag - tareDrag
    angle = np.asarray(rawdata[b"/'Untitled'/'Untitled'"])
    angle = angle[0:len(Tarm)]
    rpm = np.zeros(len(Tarm))
    t = np.arange(0, float(len(Tarm))/2000, 0.0005)
    rpm[0:len(angle)-1] = np.diff(angle)/0.0005/6
    rpm = smooth(rpm, 50)
    tsr = rpm*2*np.pi/60*0.5
    Ttare = 0.0332496307998*tsr + 0.465846267394 
    Ttrans = Ttrans + Ttare
    Tarm = Tarm + Ttare
    return t, angle, Ttrans, Tarm, drag, rpm, tsr
    
def find_t2(t, angle, t1, t2):
    angle1 = angle[2000*t1]
    angle2 = angle[2000*t2]
    N3rdRevs = np.floor((angle2-angle1)/120)
    nrevs = np.floor((angle2-angle1)/360)
    angle2 = angle1 + N3rdRevs*120
#    angle2 = angle1 + nrevs*360
    t2i = np.where(np.round(angle)==np.round(angle2))
    t2 = t[t2i]
    t2 = np.round(t2[0], decimals=2)
    return t2, nrevs
    
def calc_eta2(cp, cd):
    if cd < 0.8889:
        a = (-1+np.sqrt(1-cd))/(-2)
    elif cd >= 0.8889:  
        F = 1
        a = (18*F - 20 - 3*np.sqrt(cd*(50-36*F)+12*F*(3*F-4)))/(36*F - 50)
    eta2 = cp/((1-a)*cd)
    return a, eta2
    

def batchperf(runs="all"):
    if runs == "all":
        runs = range(1,378) # 377 runs total
    t1 = 13
    t2t = 30
    cp = np.zeros(len(runs))
    cd = np.zeros(len(runs))
    ct = np.zeros(len(runs))
    tsr = np.zeros(len(runs))
    std_tsr = np.zeros(len(runs))
    std_cp = np.zeros(len(runs))
    std_cd = np.zeros(len(runs))
    std_ct = np.zeros(len(runs))
    t2 = np.zeros(len(runs))
    eta2 = np.zeros(len(runs))
    nrevs = np.zeros(len(runs))
    a = np.zeros(len(runs))
    torque_ripple = np.zeros(len(runs))
    amp_tsr = np.zeros(len(runs))
    phase_tsr = np.zeros(len(runs))
    amp_cp = np.zeros(len(runs))
    phase_cp = np.zeros(len(runs))
    amp_cd = np.zeros(len(runs))
    phase_cd = np.zeros(len(runs))
    amp_ct = np.zeros(len(runs))
    phase_ct = np.zeros(len(runs))
    for n in range(len(runs)):
        print("Processing performance data from run", runs[n], "of", \
        str(np.max(runs))+"...")
        t, angle, Ttrans, Tarm, drag, rpm, tsr_s = loadtdms(runs[n])
        t2[n], nrevs[n] = find_t2(t, angle, t1, t2t)
        tsr[n], std_tsr[n] = calcstats(tsr_s, t1, t2[n], 2000)
        cp_s = Ttrans*tsr_s/0.5/500
        cd_s = drag/500.0
        ct_s = cp_s/tsr_s
        cp[n], std_cp[n] = calcstats(cp_s, t1, t2[n], 2000)
        cd[n], std_cd[n] = calcstats(cd_s, t1, t2[n], 2000)
        ct[n], std_ct[n] = calcstats(ct_s, t1, t2[n], 2000)
        a[n], eta2[n] = calc_eta2(cp[n], cd[n])
        torque_seg = Ttrans[2000*t1:2000*t2[n]]
        torque_ripple[n] = (np.max(torque_seg) \
                           - np.min(torque_seg))/np.mean(torque_seg)
        cp_seg = cp_s[2000*t1:2000*t2[n]]
        cd_seg = cd_s[2000*t1:2000*t2[n]]
        ct_seg = ct_s[2000*t1:2000*t2[n]]
        tsr_seg = tsr_s[2000*t1:2000*t2[n]]
        angle_seg = angle[2000*t1:2000*t2[n]]
        amp_tsr[n], phase_tsr[n] = find_amp_and_phase(angle_seg, tsr_seg)
        amp_cp[n], phase_cp[n] = find_amp_and_phase(angle_seg, cp_seg)
        amp_cd[n], phase_cd[n] = find_amp_and_phase(angle_seg, cd_seg)
        amp_ct[n], phase_ct[n] = find_amp_and_phase(angle_seg, ct_seg)
    np.save("Processed/cp", cp)
    np.save("Processed/cd", cd)
    np.save("Processed/ct", ct)
    np.save("Processed/tsr", tsr)
    np.save("Processed/std_tsr", std_tsr)
    np.save("Processed/std_cp", std_cp)
    np.save("Processed/t2", t2)
    np.save("Processed/eta2", eta2)
    np.save("Processed/nrevs", nrevs)
    np.save("Processed/a", a)
    np.save("Processed/torque_ripple", torque_ripple)
    np.save("Processed/amp_tsr", amp_tsr)
    np.save("Processed/phase_tsr", phase_tsr)
    np.save("Processed/amp_cp", amp_cp)
    np.save("Processed/phase_cp", phase_cp)
    np.save("Processed/amp_cd", amp_cd)
    np.save("Processed/phase_cd", phase_cd)
    np.save("Processed/amp_ct", amp_ct)
    np.save("Processed/phase_ct", phase_ct)
    
def find_amp_and_phase(angle, data, npeaks=3):
    amp = (np.max(data) - np.min(data))/2
    phase = angle[np.where(data == data.max())[0][0]] % (360/npeaks)
    return amp, phase

def ens_ave():
    run = 359
    t1 = 13
    t2 = 30
    t, angle, Ttrans, Tarm, drag, rpm, tsr = loadtdms(run)
    angle1 = angle[t1*2000]
    undershoot = 360-np.mod(angle1,360)
    angle1 = angle1 + undershoot
    angle2 = angle[t2*2000]
    overshoot = np.mod(angle2,360)
    angle2 = angle2 - overshoot
    nrevs = (angle2-angle1)/360
    def findIndex(angle, ta):
        i = np.where(np.round(angle, decimals=0)==ta)
        i = i[0]
        if len(i) > 1: i = i[0]
        return i
    i1 = findIndex(angle, angle1)
    i2 = findIndex(angle, angle1+360)
    npoints = i2-i1
    Tens = Ttrans[i1:i2]
    print(nrevs)
    for n in range(1,int(nrevs)):
            ta1 = angle1+n*360
            i1 = findIndex(angle, ta1)
            i2 = i1+npoints
            Tens = Tens + Ttrans[i1:i2]
    Tens = Tens/nrevs
    angleb = np.linspace(0, 360, num=npoints)
    plt.close("all")
    plt.plot(angleb, Tens, "k")
    plt.xlabel(r"$\theta$ (deg)")
    plt.ylabel(r"Torque (Nm)")
    styleplot()

def getruns(z_H, tsr):
    if z_H == 0:
        runs = range(122,167)
    if z_H == 0.25 and tsr == 1.9:
        runs = range(32, 77)
    if z_H == 0.375:
        runs = range(167,212)
    if z_H == 0.625:
        runs = range(212,257)
    if z_H == 0.125:
        runs = range(257,302)
    if z_H == 0.25 and tsr == 1.4:
        runs = range(302,347)
    if z_H == 0.5:
        runs = range(77,122)
    return runs
    
def find_run_ind(y_R, z_H, tsr):
    """Finds the run index corresponding to the inputs."""
    tp = import_testplan()
    i = np.where(np.logical_and(tp["y/R"]==y_R, 
                                tp["z/H"]==z_H,
                                tp["tsr"]==tsr))[0][0]
    return i

def batchwake():
    runs = range(1, 378)
    y_R = np.hstack([-3.,-2.75,-2.5,-2.25,-2.,-1.8,np.arange(-1.6,0.1,0.1)])
    y_R = np.hstack([y_R, -np.flipud(y_R[0:-1])])
    t1 = 13
    t2 = np.load("Processed/t2.npy")
    tsr = np.load("Processed/tsr.npy")
    meanu = np.zeros(len(runs))
    meanv = np.zeros(len(runs))
    meanw = np.zeros(len(runs))
    stdu = np.zeros(len(runs))
    stdv = np.zeros(len(runs))
    stdw = np.zeros(len(runs))
    meanuv = np.zeros(len(runs))
    meanuw = np.zeros(len(runs))
    meanvw = np.zeros(len(runs))
    meanvv = np.zeros(len(runs))
    meanww = np.zeros(len(runs))
    meanuu = np.zeros(len(runs))
    phi = np.zeros(len(runs))
    meanu2 = np.zeros(len(runs))
    vectemp = np.zeros(len(runs))
    fpeak_u = np.zeros(len(runs))
    fstrength_u = np.zeros(len(runs))
    fpeak_v = np.zeros(len(runs))
    fstrength_v = np.zeros(len(runs))
    fpeak_w = np.zeros(len(runs))
    fstrength_w = np.zeros(len(runs))
    
    for n in range(len(runs)):
        print("Processing velocity data from run", str(runs[n]) + "...")
        tv,u,v,w = loadvec(runs[n])
        phi_s = 0.5*u*(u**2 + v**2 + w**2)
        u2 = u**2
        meanu[n], stdu[n] = calcstats(u, t1, t2[n], 200)
        meanv[n], stdv[n] = calcstats(v, t1, t2[n], 200)
        meanw[n], stdw[n] = calcstats(w, t1, t2[n], 200)
        uv = (u-meanu[n])*(v-meanv[n])
        uw = (u-meanu[n])*(w-meanw[n])
        vw = (v-meanv[n])*(w-meanw[n])
        vv = (v-meanv[n])*(v-meanv[n])
        ww = (w-meanw[n])*(w-meanw[n])
        uu = (u-meanu[n])*(u-meanu[n])
        meanuv[n] = np.mean(uv[t1*200:t2[n]*200])
        meanuw[n] = np.mean(uw[t1*200:t2[n]*200])
        meanvw[n] = np.mean(vw[t1*200:t2[n]*200])
        meanvv[n] = np.mean(vv[t1*200:t2[n]*200])
        meanww[n] = np.mean(ww[t1*200:t2[n]*200])
        meanuu[n] = np.mean(uu[t1*200:t2[n]*200])
        phi[n] = np.mean(phi_s[t1*200:t2[n]*200])
        meanu2[n] = np.mean(u2[t1*200:t2[n]*200])
        vectemp[n] = loadvectemp(runs[n])
        # Spectral calculations
        u_seg = u[200*t1:200*t2[n]] - np.mean(u[200*t1:200*t2[n]])
        v_seg = v[200*t1:200*t2[n]] - np.mean(v[200*t1:200*t2[n]])
        w_seg = w[200*t1:200*t2[n]] - np.mean(w[200*t1:200*t2[n]])
        f_turbine = tsr[n]*U/R/(2*np.pi)
        # Find maximum frequency and its relative strength
        f, spec = psd(tv, u_seg, window=None)
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        fstrength_u[n] = np.max(spec)/np.var(u_seg)*(f[1] - f[0])
        fpeak_u[n] = f_max/f_turbine
        f, spec = psd(tv, v_seg, window=None)
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        fstrength_v[n] = np.max(spec)/np.var(v_seg)*(f[1] - f[0])
        fpeak_v[n] = f_max/f_turbine
        f, spec = psd(tv, w_seg, window=None)
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        fstrength_w[n] = np.max(spec)/np.var(w_seg)*(f[1] - f[0])
        fpeak_w[n] = f_max/f_turbine
    np.save("Processed/meanu", meanu)
    np.save("Processed/meanv", meanv)
    np.save("Processed/meanw", meanw)
    np.save("Processed/stdu", stdu)
    np.save("Processed/stdv", stdv)
    np.save("Processed/stdw", stdw)
    np.save("Processed/meanuv", meanuv)
    np.save("Processed/meanuw", meanuw)
    np.save("Processed/meanvw", meanvw)
    np.save("Processed/phi", phi)
    np.save("Processed/meanu2", meanu2)
    np.save("Processed/meanvv", meanvv)
    np.save("Processed/meanww", meanww)
    np.save("Processed/meanuu", meanuu)
    np.save("Processed/vectemp", vectemp)
    np.save("Processed/fpeak_u", fpeak_u)
    np.save("Processed/fstrength_u", fstrength_u)
    np.save("Processed/fpeak_v", fpeak_v)
    np.save("Processed/fstrength_v", fstrength_v)
    np.save("Processed/fpeak_w", fpeak_w)
    np.save("Processed/fstrength_w", fstrength_w)
    
def export_perf_csv(rev=0):
    """Export processed data to csv file."""
    if not os.path.isdir("Processed/csv"):
        os.mkdir("Processed/csv")
    import datetime
    datestring = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    sep = "----------------------------------------------------"
    runs = range(31)
    tsr = np.load("Processed/tsr.npy")[runs]
    tsr = np.round(tsr, decimals=2)
    cp = np.load("Processed/cp.npy")[runs]
    cp = np.round(cp, decimals=3)
    cd = np.load("Processed/cd.npy")[runs]
    cd = np.round(cd, decimals=3)
    u = np.ones(len(runs))
    vectemp = np.mean(np.load("Processed/vectemp.npy")[runs])
    vectemp = np.round(vectemp, decimals=2)
    runnumber = range(1, 32)
    metadata = ["Processed performance data from UNH-RVAT tow tank experiments",
                "Performed on March 9-10, 2013",
                "Rev" + str(rev) + " - Generated " + datestring,
                sep,
                "Turbine diameter: 1 m",
                "Turbine height: 1 m",
                "Turbine blades: 3X NACA 0020 X 0.14 m chord",
                "Turbine blade mounting: 1/2 chord, 1/2 span, zero pitch",
                "Mean ADV temperature: " + str(vectemp) + " degrees C",
                sep,
                "For more experimental details refer to:",
                'Bachant and Wosnik (2013) "Performance and Wake Measurements for a Vertical Axis Turbine at Moderate Reynolds Number"', 
                "Proceedings of ASME Fluids Engineering Division Summer Meeting 2013, Paper FED2013-16575",
                sep,
                "Columns:",
                "Run - Run number corresponding to experimental program",
                "U - Tow carriage velocity in m/s",
                "TSR - Mean turbine tip speed ratio",
                "C_P - Mean turbine power coefficient",
                "C_D - Mean turbine drag (or thrust) coefficient",
                sep]
    s = 6
    clabels = ["Run".center(s), "U".center(s), "TSR".center(s),
               "C_P".center(s), "C_D".center(s)]
    with open("Processed/csv/unh-rvat-perf-2013-03-rev"+str(rev)+".csv","wb") \
    as csvfile:
        fwriter = csv.writer(csvfile)
        for i in range(len(metadata)):
            fwriter.writerow([metadata[i]])
        fwriter.writerow(clabels)
        for run in runs:
            fwriter.writerow([str(runnumber[run]).center(s), 
                             str(u[run]).center(s), 
                             ("%.2f" %tsr[run]).center(s), 
                             ("%.3f" %cp[run]).center(s), 
                             ("%.3f" %cd[run]).center(s)])
    
def main():
    """Main function."""
#    batchperf()
#    batchwake()
#    export_perf_csv(rev=1)
    loadtdms(1)

    
if __name__ == "__main__":
    main()