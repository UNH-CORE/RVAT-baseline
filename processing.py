"""
This script imports and processes turbine performance and wake data from
March 2013 experiments with the VAT

"""
from __future__ import division, print_function
import xlrd
from pxl import timeseries, fdiff
from pxl.timeseries import *
import numpy as np
import csv
import matplotlib.pyplot as plt
import time
import matplotlib
from scipy.signal import decimate
from scipy.interpolate import interp1d
import sys
import os
import pandas as pd
import urllib
import json

try:
    import pytdms
except ImportError:
    pytdms = False

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
    return pd.read_csv("Test plan/Test plan.csv")
    
def loadvec(run):
    if not os.path.isfile("Raw/Vectrino/vec{}.dat".format(run)):
        download_vecdata(run)
    data = np.loadtxt("Raw/Vectrino/vec{}.dat".format(run))
    t = data[:,0]
    u = data[:,3]
    v = data[:,4]
    w = data[:,6]
    return t, u, v, w

def loadvectemp(run):
    if not os.path.isfile("Raw/Vectrino/vec{}.hdr".format(run)):
        download_vecdata(run)
    with open("Raw/Vectrino/vec" + str(run) + ".hdr") as f:
        temp = f.readlines()[117].split()[1]
    return float(temp)
    
def download_vecdata(run):
    """Downloads Vectrino header and data files for a given run."""
    with open("Raw/urls.json") as f:
        urls = json.load(f)
    if not os.path.isdir("Raw/Vectrino"):
        os.mkdir("Raw/Vectrino")
    print("Downloading Vectrino data file from run {}...".format(run))
    urllib.urlretrieve(urls["vec{}.dat".format(run)], 
                            filename="Raw/Vectrino/vec{}.dat".format(run))
    print("Done")
    print("Downloading Vectrino header file from run {}...".format(run))
    urllib.urlretrieve(urls["vec{}.hdr".format(run)], 
                       filename="Raw/Vectrino/vec{}.hdr".format(run))
    print("Done")
                          
def loadtdms(run):
    filename = "Raw/TDMS/run{}.tdms".format(run)
    if not os.path.isfile(filename):
        download_tdms(run)
    objects, rawdata = pytdms.read(filename)
    Ttrans = np.asarray(rawdata[b"/'Untitled'/'TorqueTrans'"])
    Tarm = np.asarray(rawdata[b"/'Untitled'/'TorqueArm'"])
    dragL = np.asarray(rawdata[b"/'Untitled'/'DragL'"])
    dragR = np.asarray(rawdata[b"/'Untitled'/'DragR'"])
    Tarm = Tarm - np.mean(Tarm[0:2000])
    taredrag = 49.2407
    drag = dragL + dragR
    drag = drag - np.mean(drag[0:2000])
    drag = drag - taredrag
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
    
def download_tdms(run):
    """This function downloads a TDMS file from figshare"""
    with open("Raw/urls.json") as f:
        urls = json.load(f)
    print("Downloading TDMS file from run {}...".format(run))
    if not os.path.isdir("Raw/TDMS"):
        os.mkdir("Raw/TDMS")
    urllib.urlretrieve(urls["run{}.tdms".format(run)], 
                       filename="Raw/TDMS/run{}.tdms".format(run))
    print("Done")
    
def find_t2(t, angle, t1, t2):
    angle1 = angle[2000*t1]
    angle2 = angle[2000*t2]
    nbladepass = np.floor((angle2-angle1)/120)
    nrevs = np.floor((angle2-angle1)/360)
    angle2 = angle1 + nbladepass*120
#    angle2 = angle1 + nrevs*360
    t2i = np.where(np.round(angle)==np.round(angle2))
    t2 = t[t2i]
    t2 = np.round(t2[0], decimals=2)
    return t2, nrevs, nbladepass
    
def calc_eta2(cp, cd):
    if cd < 0.8889:
        a = (-1+np.sqrt(1-cd))/(-2)
    elif cd >= 0.8889:  
        F = 1
        a = (18*F - 20 - 3*np.sqrt(cd*(50-36*F)+12*F*(3*F-4)))/(36*F - 50)
    eta2 = cp/((1-a)*cd)
    return a, eta2
    

def batchperf(t1=13, t2_guess=30):
    testplan = pd.read_csv("Test plan/Test plan.csv")
    try:
        df = pd.read_csv("Processed/processed.csv")
    except IOError:
        df = pd.DataFrame()
    runs = testplan["Run"]
    df["run"] = runs
    df["t1_perf"] = t1
    quantities = ["t2", "nbladepass", "tsr", "cp", "cd", "ct", "std_tsr", 
                  "std_cp", "std_cd", "std_ct", "eta2", "a", "torque_ripple",
                  "amp_tsr", "phase_tsr", "amp_cp", "phase_cp", "amp_cd",
                  "phase_cd", "amp_ct", "phase_ct"]
    for q in quantities:
        if not q in df:
            df[q] = np.zeros(len(runs))
    for n in range(len(runs)):
        print("Processing performance data from run", runs[n], "of", \
        str(np.max(runs))+"...")
        t, angle, Ttrans, Tarm, drag, rpm, tsr_s = loadtdms(runs[n])
        t2, nrevs, df.nbladepass[n] = find_t2(t, angle, t1, t2_guess)
        df.t2[n] = t2
        df.tsr[n], df.std_tsr[n] = calcstats(tsr_s, t1, t2, 2000)
        cp_s = Ttrans*tsr_s/0.5/500
        cd_s = drag/500.0
        ct_s = cp_s/tsr_s
        df.cp[n], df.std_cp[n] = calcstats(cp_s, t1, t2, 2000)
        df.cd[n], df.std_cd[n] = calcstats(cd_s, t1, t2, 2000)
        df.ct[n], df.std_ct[n] = calcstats(ct_s, t1, t2, 2000)
        df.a[n], df.eta2[n] = calc_eta2(df.cp[n], df.cd[n])
        torque_seg = Ttrans[2000*t1:2000*t2]
        df.torque_ripple[n] = (np.max(torque_seg) \
                               - np.min(torque_seg))/np.mean(torque_seg)
        cp_seg = cp_s[2000*t1:2000*t2]
        cd_seg = cd_s[2000*t1:2000*t2]
        ct_seg = ct_s[2000*t1:2000*t2]
        tsr_seg = tsr_s[2000*t1:2000*t2]
        angle_seg = angle[2000*t1:2000*t2]
        df.amp_tsr[n], df.phase_tsr[n] = find_amp_and_phase(angle_seg, tsr_seg)
        df.amp_cp[n], df.phase_cp[n] = find_amp_and_phase(angle_seg, cp_seg)
        df.amp_cd[n], df.phase_cd[n] = find_amp_and_phase(angle_seg, cd_seg)
        df.amp_ct[n], df.phase_ct[n] = find_amp_and_phase(angle_seg, ct_seg)
    # Save to CSV
    df.to_csv("Processed/processed.csv", index=False)

    
def find_amp_and_phase(angle, data, npeaks=3):
    amp = (np.max(data) - np.min(data))/2
    phase = angle[np.where(data == data.max())[0][0]] % (360/npeaks)
    return amp, phase

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
    i = np.where(np.logical_and(tp["y/R"].values==y_R, 
                                tp["z/H"].values==z_H,
                                tp["TSR"].values==tsr))[0][0]
    return i

def batchwake(t1=13):
    try:
        df = pd.read_csv("Processed/processed.csv")
        for key in df.keys():
            if "Unnamed" in key:
                del df[key]
        tsr = df["tsr"]
        t2 = df["t2"]
    except IOError:
        df = pd.DataFrame()
        t2 = np.load("Processed/t2.npy")
        tsr = np.load("Processed/tsr.npy")
    except KeyError:
        t2 = np.load("Processed/t2.npy")
        tsr = np.load("Processed/tsr.npy")
    testplan = pd.read_csv("Test plan/Test plan.csv")
    df["run"] = testplan["Run"]
    df["y/R"] = testplan["y/R"]
    df["z/H"] = testplan["z/H"]
    df["t1_vec"] = t1
    df["t2"] = t2
    runs = testplan["Run"]
    quantities = ["meanu", "meanv", "meanw", "stdu", "stdv", "stdw",
                  "meanupvp", "meanupwp", "meanvpwp", "meanupup", "meanvpvp", 
                  "meanwpwp", "meanuu", "vectemp", "fpeak_u", "fstrength_u", 
                  "fpeak_v", "fstrength_v", "fpeak_w", "fstrength_w"]
    for q in quantities:
        if not q in df:
            df[q] = np.zeros(len(runs))
    
    for n in range(len(runs)):
        print("Processing velocity data from run", str(runs[n]) + "...")
        tv,u,v,w = loadvec(runs[n])
        phi_s = 0.5*u*(u**2 + v**2 + w**2)
        uu = u**2
        meanu, stdu = calcstats(u, t1, t2[n], 200)
        meanv, stdv = calcstats(v, t1, t2[n], 200)
        meanw, stdw = calcstats(w, t1, t2[n], 200)
        df.meanu[n], df.meanv[n], df.meanw[n] = meanu, meanv, meanw
        df.stdu[n], df.stdv[n], df.stdw[n] = stdu, stdv, stdw
        upvp = (u-meanu)*(v-meanv)
        upwp = (u-meanu)*(w-meanw)
        vpwp = (v-meanv)*(w-meanw)
        vpvp = (v-meanv)*(v-meanv)
        wpwp = (w-meanw)*(w-meanw)
        upup = (u-meanu)*(u-meanu)
        df.meanupvp[n] = np.mean(upvp[t1*200:t2[n]*200])
        df.meanupwp[n] = np.mean(upwp[t1*200:t2[n]*200])
        df.meanvpwp[n] = np.mean(vpwp[t1*200:t2[n]*200])
        df.meanvpvp[n] = np.mean(vpvp[t1*200:t2[n]*200])
        df.meanwpwp[n] = np.mean(wpwp[t1*200:t2[n]*200])
        df.meanupup[n] = np.mean(upup[t1*200:t2[n]*200])
#        phi[n] = np.mean(phi_s[t1*200:t2[n]*200])
        df.meanuu[n] = np.mean(uu[t1*200:t2[n]*200])
        df.vectemp[n] = loadvectemp(runs[n])
        # Spectral calculations
        u_seg = u[200*t1:200*t2[n]] - np.mean(u[200*t1:200*t2[n]])
        v_seg = v[200*t1:200*t2[n]] - np.mean(v[200*t1:200*t2[n]])
        w_seg = w[200*t1:200*t2[n]] - np.mean(w[200*t1:200*t2[n]])
        f_turbine = tsr[n]*U/R/(2*np.pi)
        # Find maximum frequency and its relative strength
        f, spec = psd(tv, u_seg, window=None)
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        df.fstrength_u[n] = np.max(spec)/np.var(u_seg)*(f[1] - f[0])
        df.fpeak_u[n] = f_max/f_turbine
        f, spec = psd(tv, v_seg, window=None)
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        df.fstrength_v[n] = np.max(spec)/np.var(v_seg)*(f[1] - f[0])
        df.fpeak_v[n] = f_max/f_turbine
        f, spec = psd(tv, w_seg, window=None)
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        df.fstrength_w[n] = np.max(spec)/np.var(w_seg)*(f[1] - f[0])
        df.fpeak_w[n] = f_max/f_turbine
    # Save to CSV
    df.to_csv("Processed/processed.csv", index=False)
    
def convert_npy_to_csv():
    files = os.listdir("Processed")
    df = pd.DataFrame()
    for f in files:
        if (".npy") in f:
            df[f.replace(".npy", "")] = np.load("Processed/" + f)
    df.to_csv("Processed/processed.csv")
    
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
#    convert_npy_to_csv()
    
if __name__ == "__main__":
    main()