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

# Some constants
R = 0.5
H = 1.0
U = 1.0
d_shaft = 0.095
A_t = 1.0
rho = 1000.0
nu = 1e-6

def styleplot():
    font = {'family':'serif','serif':'cmr10','size':23}
    lines = {'markersize':9, 'markeredgewidth':0.9}
    legend = {'numpoints':1, 'fontsize': 'small'}
    matplotlib.rc('text', usetex = True)
    matplotlib.rc('font', **font)
    matplotlib.rc('lines', **lines)
    matplotlib.rc('legend', **legend)
    matplotlib.rc('xtick', **{'major.pad':12})
    plt.grid(True)
    plt.tight_layout()
    
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
    (objects,rawdata) = pyTDMS.read(filename)
    Ttrans = np.asarray(rawdata["/'Untitled'/'TorqueTrans'"])
    Tarm = np.asarray(rawdata["/'Untitled'/'TorqueArm'"])
    dragL = np.asarray(rawdata["/'Untitled'/'DragL'"])
    dragR = np.asarray(rawdata["/'Untitled'/'DragR'"])
    Ttrans = Ttrans - np.mean(Ttrans[0:2000])
    Tarm = Tarm - np.mean(Tarm[0:2000])
    dragL = dragL - np.mean(dragL[0:2000])
    dragR = dragR - np.mean(dragR[0:2000])
    tareDrag = 49.2407
    drag = dragL + dragR - tareDrag
    angle = np.asarray(rawdata["/'Untitled'/'Untitled'"])
    angle = angle[0:len(Tarm)]
    rpm = np.zeros(len(Tarm))
    t = np.arange(0,float(len(Tarm))/2000,0.0005)
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
    Nrevs = np.floor((angle2-angle1)/360)
    angle2 = angle1 + N3rdRevs*120
#    angle2 = angle1 + Nrevs*360
    t2i = np.where(np.round(angle)==np.round(angle2))
    t2 = t[t2i]
    t2 = np.round(t2[0], decimals=2)
    return t2, Nrevs
    
def calc_eta2(cp, cd):
    if cd < 0.8889:
        a = (-1+np.sqrt(1-cd))/(-2)
    if cd >= 0.8889:  
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
    tsr = np.zeros(len(runs))
    std_tsr = np.zeros(len(runs))
    std_cp = np.zeros(len(runs))
    std_cd = np.zeros(len(runs))
    t2 = np.zeros(len(runs))
    eta2 = np.zeros(len(runs))
    Nrevs = np.zeros(len(runs))
    a = np.zeros(len(runs))
    torque_ripple = np.zeros(len(runs))
    for n in range(len(runs)):
        print("Processing performance data from run", runs[n], "of", \
        str(np.max(runs))+"...")
        t, angle, Ttrans, Tarm, drag, rpm, tsr_s = loadtdms(runs[n])
        t2[n], Nrevs[n] = find_t2(t, angle, t1, t2t)
        tsr[n], std_tsr[n] = calcstats(tsr_s, t1, t2[n], 2000)
        cp_s = Ttrans*tsr_s/0.5/500
        cp[n], std_cp[n] = calcstats(cp_s, t1, t2[n], 2000)
        cd[n], std_cd[n] = calcstats(drag/500, t1, t2[n], 2000)
        a[n], eta2[n] = calc_eta2(cp[n], cd[n])
        torque_seg = Ttrans[2000*t1:2000*t2[n]]
        torque_ripple[n] = (np.max(torque_seg) \
                           - np.min(torque_seg))/np.mean(torque_seg)
    np.save('Processed/cp', cp)
    np.save('Processed/cd', cd)
    np.save('Processed/tsr', tsr)
    np.save('Processed/std_tsr', std_tsr)
    np.save('Processed/std_cp', std_cp)
    np.save('Processed/t2', t2)
    np.save('Processed/eta2', eta2)
    np.save('Processed/Nrevs', Nrevs)
    np.save('Processed/a', a)
    np.save('Processed/torque_ripple', torque_ripple)

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
    Nrevs = (angle2-angle1)/360
    def findIndex(angle, ta):
        i = np.where(np.round(angle, decimals=0)==ta)
        i = i[0]
        if len(i) > 1: i = i[0]
        return i
    i1 = findIndex(angle, angle1)
    i2 = findIndex(angle, angle1+360)
    npoints = i2-i1
    Tens = Ttrans[i1:i2]
    print(Nrevs)
    for n in range(1,int(Nrevs)):
            ta1 = angle1+n*360
            i1 = findIndex(angle, ta1)
            i2 = i1+npoints
            Tens = Tens + Ttrans[i1:i2]
    Tens = Tens/Nrevs
    angleb = np.linspace(0, 360, num=npoints)
    plt.close('all')
    plt.plot(angleb, Tens, 'k')
    plt.xlabel(r'$\theta$ (deg)')
    plt.ylabel(r'Torque (Nm)')
    styleplot()

def plotsinglerun(run, perf=True, wake=False, autocorr=False, save=False):
    t1 = 13
    t2 = 30
    t2t = 30
    t, angle, Ttrans, Tarm, drag, rpm, tsr = loadtdms(run)
    t2, Nrevs = find_t2(t, angle, t1, t2t)
    meantsr, std_tsr = calcstats(tsr, t1, t2, 2000)
    omega = meantsr*U/R
    blade_period = 2*np.pi/omega/3
    vecloaded = False
    if perf:
        cp = (Ttrans+0.5)*tsr/0.5/500.0
        cd_ts = drag/500.0
        max_cd = np.max(cd_ts[t1*2000:t2*2000])
        min_cd = np.min(cd_ts[t1*2000:t2*2000])
        max_torque = np.max(Ttrans[t1*2000:t2*2000])
        min_torque = np.min(Ttrans[t1*2000:t2*2000])
        cd, std_cd = calcstats(cd_ts, t1, t2, 2000)
        cp, std_cp = calcstats(cp, t1, t2, 2000)
        meanT, stdT = calcstats(Ttrans, t1, t2, 2000)
        meanrpm, std_rpm = calcstats(rpm, t1, t2, 2000)
#        plt.plot(t, Ttrans, 'k')
        plt.hold(True)
        plt.plot(t, Ttrans, 'k')
        plt.xlabel(r'$t$(s)')
#        plt.ylabel(r'Torque (Nm)')
#        plt.vlines([t1, t2],np.min(Ttrans),np.max(Ttrans),
#                   color='r',linestyles='dashed')
#        plt.vlines([t2t],np.min(Ttrans),np.max(Ttrans),
#                   color='k',linestyles='dashed')
        styleplot()
        print('TSR =', meantsr, '; C_P =', cp)
        print("Number of revolutions:", Nrevs)
        print("(1/2)(max_cd-min_cd)/cd:", np.abs(max_cd-min_cd)*0.5/cd)
        print("(1/2)(Tmax-Tmin)/Tmean:", np.abs(max_torque-min_torque)*0.5/meanT)
    if wake:
        tv, u, v, w = loadvec(run)
        vecloaded = True
        angle = decimate(angle, 10)
        meanu, stdu = calcstats(u, t1, t2, 200)
        plt.figure()
        plt.plot(tv, u, 'k')
        plt.xlabel('$t$ (s)')
        plt.ylabel('$u$ (m/s)')
        plt.vlines([t1, t2],np.min(u),np.max(u),
                   color='r',linestyles='dashed')
        styleplot()
        plt.figure()
        plt.plot(angle, u[:len(angle)])
    if autocorr:
        if not vecloaded:
            tv, u, v, w = loadvec(run)
        u = u[t1*200:t2*200]
        t = tv[t1*200:t2*200]
        # Compute autocorrelation
        tau, rho = timeseries.autocorr(u, t, 0, 6.0)
        print("Blade passage period =", blade_period, "s")
        # Compute integral timescale for velocity
        i = np.where(np.round(rho, decimals=2)==0)[0][0]
        int_time = tau[i]
        print("Integral timescale =", int_time)
        plt.figure()
        plt.plot(tau, rho)
        plt.vlines([blade_period], rho.min(), 1,
                   color='k',linestyles='dashed')
        plt.vlines([blade_period*3], rho.min(), 1,
                   color='k',linestyles='dashed')
        plt.xlabel("Lag (s)")
        plt.ylabel("Autocorrelation coefficient")
        styleplot()
    plt.show()

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

def plotvelspec(y_R=0, z_H=0, tsr=1.9, show=False):
    """Plots the velocity spectrum for a single run."""
    # Find index for the desired parameters
    i = find_run_ind(y_R, z_H, tsr)
    print("Plotting spectra from run", i+1)
    t1 = 13
    t2 = np.load('Processed/t2.npy')[i]
    t, u, v, w = loadvec(i+1) # Run name is index + 1
    v_seg = v[200*t1:200*t2] - np.mean(v[200*t1:200*t2])
    f, spec = psd(t, v_seg, window="Hanning")
    f_turbine = tsr*U/R/(2*np.pi)
    f_blade = f_turbine*3
    # Find maximum frequency and its relative strength
    f_max = f[np.where(spec==np.max(spec))[0][0]]
    strength = np.max(spec)/np.var(v_seg)
    if show:
        print("Strongest frequency f/f_turbine:", f_max/f_turbine)
        print("Relative strength:", strength)
        # Calculate shaft shedding frequency
        St = 0.19 # Approximate for Re_d = 1e5
        f_cyl = St*U/d_shaft
        plt.figure()
        plt.loglog(f/f_turbine, spec, 'k')
        plt.xlim((0, 50))
        plt.xlabel(r"$f/f_{\mathrm{turbine}}$")
        plt.ylabel(r"Power spectral density")
        # Should the spectrum be normalized?
        plot_vertical_lines([1, 3])
        plt.ylim((1e-6, 1e-1))
        f_line = np.linspace(10,40)
        spec_line = f_line**(-5./3)*0.05
        plt.hold(True)
        plt.loglog(f_line, spec_line)
        styleplot()
        plt.grid()
        plt.show()
    
def plot_vertical_lines(x):
    ymin = plt.axis()[2]
    ymax = plt.axis()[3]*100
    plt.vlines(x, ymin, ymax,
               color='gray', linestyles='dashed')
    plt.ylim((ymin, ymax))
    
def plotvelhist(run):
    """Plots the velocity histogram for a given run."""
    i = run - 1 # Run indexing starts from 1!
    t1 = 13
    t2 = np.load("Processed/t2.npy")[i]
    t, u, v, w = loadvec(run)
    u = u[t1*200:t2*200]
    plt.figure()
    plt.hist(u, bins=50, histtype="step", color="k", normed=True)
    plt.xlabel(r"$u/U_\infty$")
    plt.ylabel("Samples (normalized)")
    styleplot()
    plt.grid(False)
    plt.show()

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
    fpeak = np.zeros(len(runs))
    fstrength = np.zeros(len(runs))
    
    for n in range(len(runs)):
        print("Processing velocity data from run", runs[n])
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
        v_seg = v[200*t1:200*t2[n]] - np.mean(v[200*t1:200*t2[n]])
        f, spec = psd(tv, v_seg, window="Hanning")
        f_turbine = tsr[n]*U/R/(2*np.pi)
        # Find maximum frequency and its relative strength
        f_max = f[np.where(spec==np.max(spec))[0][0]]
        fstrength[n] = np.max(spec)/np.var(v_seg)*(f[1] - f[0])
        fpeak[n] = f_max/f_turbine
    np.save('Processed/meanu', meanu)
    np.save('Processed/meanv', meanv)
    np.save('Processed/meanw', meanw)
    np.save('Processed/stdu', stdu)
    np.save('Processed/stdv', stdv)
    np.save('Processed/stdw', stdw)
    np.save('Processed/meanuv', meanuv)
    np.save('Processed/meanuw', meanuw)
    np.save('Processed/meanvw', meanvw)
    np.save('Processed/phi', phi)
    np.save('Processed/meanu2', meanu2)
    np.save('Processed/meanvv', meanvv)
    np.save('Processed/meanww', meanww)
    np.save('Processed/meanuu', meanuu)
    np.save('Processed/vectemp', vectemp)
    np.save('Processed/fpeak', fpeak)
    np.save('Processed/fstrength', fstrength)

def plotperf(save=False, savepath="", savetype=".pdf"):
    i = range(31)
    cp = np.load('Processed/cp.npy')
    cd = np.load('Processed/cd.npy')
    tsr = np.load('Processed/tsr.npy')
    eta2 = np.load('Processed/eta2.npy')
    # Exergy efficiency plot
    plt.figure()
    plt.plot(tsr[i], eta2[i], '-ok', markerfacecolor = 'none')
    plt.xlabel(r'$\lambda$', labelpad=20)
    plt.ylabel(r'$\eta_{II}$')
    styleplot()
    if save:
        plt.savefig(savepath+'eta2'+savetype)
    # Power coefficient plot
    plt.figure()
    plt.plot(tsr[i], cp[i], '-ok', markerfacecolor = 'none')
    plt.xlabel(r'$\lambda$', labelpad=20)
    plt.ylabel(r'$C_P$')
    styleplot()
    if save:
        plt.savefig(savepath+'cpvstsr'+savetype)
    # Drag coefficient plot
    plt.figure()
    plt.plot(tsr[i], cd[i], '-ok', markerfacecolor = 'none')
    plt.xlabel(r'$\lambda$', labelpad=20)
    plt.ylabel(r'$C_D$')
    plt.ylim(0, 1.2001)
    styleplot()
    if save:
        plt.savefig(savepath+'cdvstsr'+savetype)
    # Torque coefficient plot
    ct = cp/tsr
    plt.figure()
    plt.plot(tsr[i], ct[i], '-ok', markerfacecolor = 'none')
    plt.xlabel(r'$\lambda$', labelpad=20)
    plt.ylabel(r'$C_T$')
    styleplot()
    if save:
        plt.savefig(savepath+'ctvstsr'+savetype)
    
def plotwake(plotlist, save=False, savepath=None, savetype=".pdf"):
    z_H = np.arange(0, 0.75, 0.125)
    y_R = np.hstack([-3.,-2.75,-2.5,-2.25,-2.,-1.8,np.arange(-1.6,0.1,0.1)])
    y_R = np.hstack([y_R, -np.flipud(y_R[0:-1])])
    y_R = np.round(y_R, decimals=4)
    # Load processed data
    meanu = np.load('Processed/meanu.npy')
    meanv = np.load('Processed/meanv.npy')
    meanw = np.load('Processed/meanw.npy')
    stdu = np.load('Processed/stdu.npy')
    stdv = np.load('Processed/stdv.npy')
    stdw = np.load('Processed/stdw.npy')
    meanuv = np.load('Processed/meanuv.npy')
    meanuw = np.load('Processed/meanuw.npy')
    meanvw = np.load('Processed/meanvw.npy')
    meanvv = np.load('Processed/meanvv.npy')
    meanww = np.load('Processed/meanww.npy')
    meanuu = np.load('Processed/meanuu.npy')
    phi = np.load('Processed/phi.npy')
    meanu2 = np.load('Processed/meanu2.npy')
    fpeak = np.load('Processed/fpeak.npy')
    fstrength = np.load('Processed/fstrength.npy')
    k = 0.5*(stdu**2 + stdv**2 + stdw**2)
    meank = 0.5*(meanu**2 + meanv**2 + meanw**2)
    kbar = meank + k
    # Create empty 2D arrays for contour plots, etc.
    meanu_a = np.zeros((len(z_H), len(y_R)))
    meanv_a = np.zeros((len(z_H), len(y_R)))
    meanw_a = np.zeros((len(z_H), len(y_R)))
    stdu_a = np.zeros((len(z_H), len(y_R)))
    stdv_a = np.zeros((len(z_H), len(y_R)))
    stdw_a = np.zeros((len(z_H), len(y_R)))
    meanuv_a = np.zeros((len(z_H), len(y_R)))
    meanuw_a = np.zeros((len(z_H), len(y_R)))
    meanvw_a = np.zeros((len(z_H), len(y_R)))
    meanvv_a = np.zeros((len(z_H), len(y_R)))
    meanww_a = np.zeros((len(z_H), len(y_R)))
    meanuu_a = np.zeros((len(z_H), len(y_R)))
    phi_a = np.zeros((len(z_H), len(y_R)))
    k_a = np.zeros((len(z_H), len(y_R)))
    meank_a = np.zeros((len(z_H), len(y_R)))
    meanu2_a = np.zeros((len(z_H), len(y_R)))
    fpeak_a = np.zeros((len(z_H), len(y_R)))
    fstrength_a = np.zeros((len(z_H), len(y_R)))
    kbar_a = np.zeros((len(z_H), len(y_R)))
    # Populate 2D arrays for velocity fields
    for n in range(len(z_H)):
        runs = getruns(z_H[n], tsr=1.9)
        ind = [run - 1 for run in runs]
        meanu_a[n,:] = meanu[ind]
        meanv_a[n,:] = meanv[ind]
        meanw_a[n,:] = meanw[ind]
        stdu_a[n,:] = stdu[ind]
        stdv_a[n,:] = stdv[ind]
        stdw_a[n,:] = stdw[ind]
        meanuv_a[n,:] = meanuv[ind]
        meanuw_a[n,:] = meanuw[ind]
        meanvw_a[n,:] = meanvw[ind]
        meanvv_a[n,:] = meanvv[ind]
        meanww_a[n,:] = meanww[ind]
        meanuu_a[n,:] = meanuu[ind]
        phi_a[n,:] = phi[ind]
        k_a[n,:] = k[ind]
        meank_a[n,:] = meank[ind]
        meanu2_a[n,:] = meanu2[ind]
        kbar_a[n,:] = kbar[ind]
        fpeak_a[n,:] = fpeak[ind]
        fstrength_a[n,:] = fstrength[ind]
    def turb_lines():
        plt.hlines(0.5, -1, 1, linestyles='solid', linewidth=2)
        plt.vlines(-1, 0, 0.5, linestyles='solid', linewidth=2)
        plt.vlines(1, 0, 0.5, linestyles='solid', linewidth=2)
    def calc_meankturbtrans():
        z = H*z_H
        y = R*y_R
        ddy_uvU = np.zeros(np.shape(meanu_a))
        ddz_uwU = np.zeros(np.shape(meanu_a))
        ddy_vvV = np.zeros(np.shape(meanu_a))
        ddz_vwV = np.zeros(np.shape(meanu_a))
        ddy_vwW = np.zeros(np.shape(meanu_a))
        ddz_wwW = np.zeros(np.shape(meanu_a))
        for n in range(len(z)):
            ddy_uvU[n,:] = fdiff.second_order_diff((meanuv_a*meanu_a)[n,:], y)
            ddy_vvV[n,:] = fdiff.second_order_diff((meanvv_a*meanv_a)[n,:], y)
            ddy_vwW[n,:] = fdiff.second_order_diff((meanvw_a*meanw_a)[n,:], y)
        for n in range(len(y)):
            ddz_uwU[:,n] = fdiff.second_order_diff((meanuw_a*meanu_a)[:,n], z)
            ddz_vwV[:,n] = fdiff.second_order_diff((meanvw_a*meanv_a)[:,n], z)
            ddz_wwW[:,n] = fdiff.second_order_diff((meanww_a*meanw_a)[:,n], z)
        tt = -0.5*(ddy_uvU + ddz_uwU + ddy_vvV + ddz_vwV + ddy_vwW + ddz_wwW)
        tty = -0.5*(ddy_uvU + ddy_vvV + ddy_vwW) # Only ddy terms
        ttz = -0.5*(ddz_uwU + ddz_vwV + ddz_wwW) # Only ddz terms
        return tt, tty, ttz
    def calc_kprod_meandiss():
        z = H*z_H
        y = R*y_R
        dUdy = np.zeros(np.shape(meanu_a))
        dUdz = np.zeros(np.shape(meanu_a))
        dVdy = np.zeros(np.shape(meanu_a))
        dVdz = np.zeros(np.shape(meanu_a))
        dWdy = np.zeros(np.shape(meanu_a))
        dWdz = np.zeros(np.shape(meanu_a))
        for n in xrange(len(z)):
            dUdy[n,:] = fdiff.second_order_diff(meanu_a[n,:], y)
            dVdy[n,:] = fdiff.second_order_diff(meanv_a[n,:], y)
            dWdy[n,:] = fdiff.second_order_diff(meanw_a[n,:], y)
        for n in xrange(len(y)):
            dUdz[:,n] = fdiff.second_order_diff(meanu_a[:,n], z)
            dVdz[:,n] = fdiff.second_order_diff(meanv_a[:,n], z)
            dWdz[:,n] = fdiff.second_order_diff(meanw_a[:,n], z)
        kprod = meanuv_a*dUdy + meanuw_a*dUdz + meanvw_a*dVdz + meanvw_a*dWdy\
                + meanvv_a*dVdy + meanww_a*dWdz
        meandiss = -2.0*nu*(dUdy**2 + dUdz**2 + dVdy**2 + dVdz**2 + dWdy**2 + dWdz**2)
        return kprod, meandiss
    def calc_meankgrad():
        z = H*z_H
        y = R*y_R
        dKdy = np.zeros(np.shape(meanu_a))
        dKdz = np.zeros(np.shape(meanu_a))
        for n in xrange(len(z)):
            dKdy[n,:] = fdiff.second_order_diff(meank_a[n,:], y)
        for n in xrange(len(y)):
            dKdz[:,n] = fdiff.second_order_diff(meank_a[:,n], z)
        return dKdy, dKdz
    if "meanucont" in plotlist or "all" in plotlist:
        # Plot contours of mean streamwise velocity
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, meanu_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.3)
        cb.set_label(r'$\overline{u}/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'meanucont'+savetype)
    if "v-wquiver" in plotlist or "all" in plotlist:
        # Make quiver plot of v and w velocities
        plt.figure(figsize=(10,5))
        Q = plt.quiver(y_R, z_H,meanv_a, meanw_a, angles='xy')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.2, 0.1, r'$0.1$ m/s',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.tight_layout()
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=2)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'v-wquiver'+savetype)
    if 'meanu_2tsrs' in plotlist:
        # Plot mean velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanu[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanu[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$\overline{u}/U_\infty$')
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+'meanu_2tsrs'+savetype)
    if 'stdu_2tsrs' in plotlist:
        # Plot stdu velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdu[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdu[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$\sigma_u/U_\infty$')
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+'stdu_2tsrs'+savetype)
    if 'uw_2tsrs' in plotlist:
        # Plot uw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuw[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuw[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r"$\overline{u'w'}$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+'uw_2tsrs'+savetype)
    if 'meanuvstsr' in plotlist:
        # Plot mean velocity components vs TSR
        tsr = np.load('Processed/tsr.npy')
        runs = range(1,32)
        ind = [run-1 for run in runs]
        plt.figure()
        plt.plot(tsr[ind], meanu[ind], '-ok', markerfacecolor='none')
        plt.xlabel(r'$\lambda$')
        plt.ylabel(r'$\overline{u}/U_\infty$')
        plt.hold(True)
        plt.plot(tsr[ind], meanv[ind], '-sk', markerfacecolor='none')
        plt.plot(tsr[ind], meanw[ind], '-^k', markerfacecolor='none')
        runs = range(347,378)
        ind = [run-1 for run in runs]
        plt.plot(tsr[ind], meanu[ind], '-ok', markerfacecolor='k',
             label=r'$\overline{u}$')
        plt.plot(tsr[ind], meanv[ind], '-sk', markerfacecolor='k',
             label=r'$\overline{v}$')
        plt.plot(tsr[ind], meanw[ind], '-^k', markerfacecolor='k',
             label=r'$\overline{w}$')
        plt.legend(ncol=3)
        styleplot()
        if save:
            plt.savefig(savepath+'meanuvstsr'+savetype)
    if 'stducont' in plotlist:
        # Plot contours of streamwise turbulence intensity
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, stdu_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.3)
        cb2.set_label(r'$\sigma_u/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'stducont'+savetype)
    if 'uvcont' in plotlist:
        # Plot contours of uv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanuv_a, 15, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                           orientation='horizontal', pad=0.26)
        cb2.set_label(r"$\overline{u'v'}/U_\infty^2$")
        cb2.set_ticks(np.arange(-0.02, 0.025, 0.005), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'uvcont'+savetype)
    if 'meanw_2tsrs' in plotlist:
        # Plot mean vertical velocity profiles at two TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanw[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanw[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$\overline{w}/U_\infty$')
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+'meanw_2tsrs'+savetype)
    if 'meanv_2tsrs' in plotlist:
        # Plot mean cross stream velocity profiles at two TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanv[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanv[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$\overline{v}/U_\infty$')
        plt.ylim(-0.1501, 0.1001)
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+'meanv_2tsrs'+savetype)
    if 'stdv_2tsrs' in plotlist:
        # Plot stdv velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdv[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdv[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$\sigma_v/U_\infty$')
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+'stdv_2tsrs'+savetype)
    if 'stdw_2tsrs' in plotlist:
        # Plot stdw velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdw[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdw[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$\sigma_w/U_\infty$')
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+'stdw_2tsrs'+savetype)
    if 'uv_2tsrs' in plotlist:
        # Plot uv Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuv[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuv[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r"$\overline{u'v'}$")
        plt.legend()
        styleplot()
        if save:
            plt.savefig(savepath+'uv_2tsrs'+savetype)
    if "kcont" in plotlist or "all" in plotlist:
        # Plot contours of k
        plt.figure(figsize=(10,5))
        csphi = plt.contourf(y_R, z_H, k_a/(0.5*1.0**2), 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cbphi = plt.colorbar(csphi, shrink=1, extend='both', 
                             orientation='horizontal', pad=0.26)
        cbphi.set_label(r'$k/\frac{1}{2}U_\infty^2$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'kcont'+savetype)
    if 'meankcont' in plotlist:
        # Plot contours of k
        plt.figure(figsize=(10,5))
        csphi = plt.contourf(y_R, z_H, meank_a/(0.5*1**2), 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cbphi = plt.colorbar(csphi, shrink=1, extend='both', 
                             orientation='horizontal', pad=0.26)
        cbphi.set_label(r'$K/\frac{1}{2}U_\infty^2$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'meankcont'+savetype)
    if 'meanvcont' in plotlist:
        # Plot contours of meanv
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, meanv_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cbmv.set_label(r'$\overline{v}/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'meanvcont'+savetype)
    if 'stdvcont' in plotlist:
        # Plot contours of stdv
        plt.figure(figsize=(10,5))
        cstdv = plt.contourf(y_R, z_H, stdv_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cbstdv = plt.colorbar(cstdv, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cbstdv.set_label(r'$\sigma_v/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'stdvcont'+savetype)
    if 'meanwcont' in plotlist:
        # Plot contours of meanw
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, meanw_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cbmv.set_label(r'$\overline{w}/U_{\infty}$')
        turb_lines()
        if save:
            plt.savefig(savepath+'meanwcont'+savetype)
    if 'stdwcont' in plotlist:
        # Plot contours of stdw
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, stdw_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cbmv.set_label(r'$\sigma_w/U_{\infty}$')
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'stdwcont'+savetype)
    if 'vw_2tsrs' in plotlist:
        # Plot vw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvw[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvw[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r"$\overline{v'w'}$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+'vw_2tsrs'+savetype)
    if 'vwcont' in plotlist:
        # Plot contours of vw Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanvw_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                           orientation='horizontal', pad=0.3)
        cb2.set_label(r"$\overline{v'w'}/U_\infty^2$")
        cb2.set_ticks(np.linspace(-.008,.006,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'vwcont'+savetype)
    if 'uwcont' in plotlist:
        # Plot contours of vw Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanuw_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                           orientation='horizontal', pad=0.26)
        cb2.set_label(r"$\overline{u'w'}/U_\infty^2$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'uwcont'+savetype)
    if 'vvcont' in plotlist:
        # Plot contours of vv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanvv_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cb2.set_label(r"$\overline{v'v'}$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'vvcont'+savetype)
    if 'wwcont' in plotlist:
        # Plot contours of vv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanww_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cb2.set_label(r"$\overline{w'w'}$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'wwcont'+savetype)
    if 'uucont' in plotlist:
        # Plot contours of uu Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanuu_a, 20)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                           orientation='horizontal', pad=0.2)
        cb2.set_label(r"$\overline{u'u'}$")
        cb2.set_ticks(np.linspace(0,.108,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'uucont'+savetype)
    if 'vv_2tsrs' in plotlist:
        # Plot vw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvv[ind], '-ok', markerfacecolor='none', 
                 label=r'$\lambda = 1.9$')
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvv[ind], '--^k', markerfacecolor='none',
                 label=r'$\lambda=1.4$')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r"$\overline{v'v'}$")
        plt.legend()
        styleplot()
        if save:
            plt.savefig(savepath+'vv_2tsrs'+savetype)
    if "fpeak" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fpeak_a, cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                           orientation='horizontal', pad=0.26)
        cb2.set_label(r"$f_{\mathrm{peak}}/f_{\mathrm{turbine}}$")
        cb2.set_ticks(np.linspace(0,10,11), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'fpeak'+savetype)
    if "fstrength" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fstrength_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend='both', 
                           orientation='horizontal', pad=0.26)
        cb2.set_label(r"$\Phi_{\max}/\sigma^2_v$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'fstrength'+savetype)
    # Plot estimate for production of turbulence kinetic energy
    if "kprod" in plotlist:
        calc_meanvelgrad()
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, kprod, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.26)
        cb.set_label(r"$-\overline{u_i' u_j'}\frac{\partial U_i}{\partial x_j}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'kprod'+savetype)
    if 'meankadv' in plotlist:
        calc_meankgrad()
        # Make quiver plot of K advection
        plt.figure(figsize=(10,5))
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='r',
                   linewidth=3)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=3)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='r',
                   linewidth=3)
        Q = plt.quiver(y_R, z_H, meanv_a/meanu_a*dKdy, meanw_a/meanu_a*dKdz, 
                       scale=4, angles='xy')
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.1, 0.2, r'$0.2 \mathrm{\, m/s^2}$',
                      labelpos='E',
                      coordinates='figure',
                      fontproperties={'size': 'small'})
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.grid(False)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'meankadv'+savetype)
    if "Kturbtrans" in plotlist:
        calc_meankturbtrans()
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, tty, 20, cmap=plt.cm.coolwarm,
                          levels=np.linspace(-0.08, 0.08, 21))
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.26)
        cb.set_label(r"$-\frac{1}{2}\frac{\partial}{\partial x_j}\overline{u_i' u_j'} U_i$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+'Kturbtrans'+savetype)
    if "meancomboquiv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,6))
        # Add contours of mean velocity
        cs = plt.contourf(y_R, z_H, meanu_a, 20, cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.2)
        cb.set_label(r'$U/U_{\infty}$')
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(y_R, z_H, meanv_a, meanw_a, angles='xy', width=0.0022)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.3, 0.1, r'$0.1 U_\infty$',
                   labelpos='E',
                   coordinates='figure',
                   fontproperties={'size': 'small'})
        plt.hlines(0.5, -1, 1, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.vlines(-1, -0.2, 0.5, linestyles='solid', colors='gray',
                   linewidth=3)
        plt.vlines(1, -0.2, 0.5, linestyles='solid', colors='gray',
                   linewidth=3)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"meancomboquiv"+savetype)
    if 'xvorticity' in plotlist:
        z = 1.0*z_H
        y = R*y_R
        dWdy = np.zeros(np.shape(meanu_a))
        dVdz = np.zeros(np.shape(meanu_a))
        for n in xrange(len(z)):
            dWdy[n,:] = fdiff.second_order_diff(meanw_a[n,:], y)
        for n in xrange(len(y)):
            dVdz[:,n] = fdiff.second_order_diff(meanv_a[:,n], z)
        # Make quiver plot of K advection
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, dWdy-dVdz, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r'$y/R$')
        plt.ylabel(r'$z/H$')
        cb = plt.colorbar(cs, shrink=1, extend='both', 
                          orientation='horizontal', pad=0.26)
        cb.set_label(r"$\Omega_x$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+'xvorticity'+savetype)
    if "Kbargraphs" in plotlist:
        """Make a bar graph of terms contributing to dK/dx:
          * Cross-stream advection
          * Vertical advection
          * Transport by turbulent fluctuations
          * Production of TKE
          * Mean dissipation
        """
        tt, tty, ttz = calc_meankturbtrans()
        kprod, meandiss = calc_kprod_meandiss()
        dKdy, dKdz = calc_meankgrad()
        plt.figure(figsize=(13,12))
        names = [r"$y$-adv.", r"$z$-adv.", 
                 r"$y$-turb.", 
                 r"$z$-turb.",
                 r"$k$-prod.", "$U$-diss."]
        locs = [(-1, 0.5, 0, 0, "(a)"),
                (-1, 0.25, 1, 0, "(b)"),
                (-1, 0.0, 2, 0, "(c)"),
                (0.0, 0.0, 2, 1, "(d)"),
                (0.0, 0.25, 1, 1, "(e)"),
                (0.0, 0.5, 0, 1, "(f)"),
                (1, 0.0, 2, 2, "(g)"),
                (1, 0.25, 1, 2, "(h)"),
                (1, 0.5, 0, 2, "(i)")]
        for yloc, zloc, ig0, ig1, letter in locs:
            i1 = np.where(z_H==zloc)[0]
            i2 = np.where(y_R==yloc)[0]
            quantities = [-meanv_a[i1,i2]/meanu_a[i1,i2]*dKdy[i1,i2], 
                          -meanw_a[i1,i2]/meanu_a[i1,i2]*dKdz[i1,i2],
                          tty[i1,i2]/meanu_a[i1,i2],
                          ttz[i1,i2]/meanu_a[i1,i2],
                          kprod[i1,i2]/meanu_a[i1,i2],
                          meandiss[i1,i2]/meanu_a[i1,i2]]
            ax = plt.subplot2grid((3,3), (ig0, ig1))
            ax.bar(range(len(names)), quantities, width=0.5)
            ax.set_xticks(np.arange(len(names))+0.25)
            ax.set_xticklabels(names, fontsize=10)
            plt.hlines(0, 0, len(names), color="gray")
            plt.title(r"$y/R =" + str(yloc) + "$; $z/H =" + str(zloc) + "$",
                      fontsize=16)
            plt.ylim((-0.2, 0.2))
        styleplot()
        plt.grid(False)
    if "Kbargraph" in plotlist:
        """Make a bar graph of terms contributing to dK/dx:
          * Cross-stream advection
          * Vertical advection
          * Transport by turbulent fluctuations
          * Production of TKE
          * Mean dissipation
        """
        tt, tty, ttz = calc_meankturbtrans()
        kprod, meandiss = calc_kprod_meandiss()
        dKdy, dKdz = calc_meankgrad()
        plt.figure(figsize=(10,5))
        names = [r"$y$-adv.", r"$z$-adv.", 
                 r"$y$-turb.", 
                 r"$z$-turb.",
                 r"$k$-prod.", "Mean diss."]
        quantities = [average_over_area(-meanv_a/meanu_a*dKdy, y_R, z_H), 
                      average_over_area(-meanw_a/meanu_a*dKdz, y_R, z_H),
                      average_over_area(tty/meanu_a, y_R, z_H),
                      average_over_area(ttz/meanu_a, y_R, z_H),
                      average_over_area(kprod/meanu_a, y_R, z_H),
                      average_over_area(meandiss/meanu_a, y_R, z_H)]
        ax = plt.gca()
        ax.bar(range(len(names)), quantities, color="k", width=0.5)
        ax.set_xticks(np.arange(len(names))+0.25)
        ax.set_xticklabels(names)
        plt.hlines(0, 0, len(names), color="gray")
        plt.ylabel(r"$\frac{K\mathrm{-transport}}{U}$ $(\mathrm{m}/\mathrm{s}^2)$")
#        ax.annotate(r"$\mathrm{Total} = " \
#                    + str(np.round(np.sum(quantities), decimals=4)) + "$", 
#                    xy=(0, 0), xytext=(0.75, 0.82), 
#                    xycoords="figure fraction", fontsize=18)
        styleplot()
        plt.grid(False)
        if save:
            plt.savefig(savepath+'Kbargraph'+savetype)
    plt.show()
    # Look at exergy efficiency -- probably wrong
    # Calculate spatial average <> of velocities and velocities squared
    Ubr_1 = 1.0 
    Ubr_2 = np.trapz(np.trapz(meanu_a, y_R*R, axis=1), dx=0.125)/(3.0*0.625)
    U2br_1 = 1.0**2
    U2br_2 = np.trapz(np.trapz(meanu2_a, y_R*R, axis=1), dx=0.125)/(3.0*0.625)
    kbarbr_1 = 0.5*1.0**2
    kbarbr_2 = np.trapz(np.trapz(kbar_a, y_R*R, axis=1), dx=0.125)/(3.0*0.625)
    phibr_1 = 0.5**1.0*1.0**2
    phibr_2 = np.trapz(np.trapz(phi_a, y_R*R, axis=1), dx=0.125)/(3.0*0.625)
    A_2 = 3*0.625*2
    A_1 = A_2*Ubr_2/Ubr_1 # Solve for A1 using continuity
    cd_meas = 0.964 # .964
    fd_meas = 0.5*rho*cd_meas*A_t*1.0**2
    p_1 = 0.0 # Gage pressure 1 in Pascals
    p_2 = -(fd_meas - p_1*A_1 - rho*(A_1*U2br_1 - A_2*U2br_2))/A_2
    power_d = rho*(Ubr_1*A_1*(kbarbr_1 + p_1/rho) - Ubr_2*A_2*(kbarbr_2 + p_2/rho))
    power_d = rho*phibr_1*A_1 + Ubr_1*p_1*A_1 - rho*phibr_2*A_2 - Ubr_2*p_2*A_2
    power_d_k = rho*phibr_1*A_1 - rho*phibr_2*A_2 
    power_d_p = Ubr_1*p_1*A_1 - Ubr_2*p_2*A_2
    ptot1 = 0.5*U**2 + 0.0
    ptot2 = ptot1 - fd_meas/rho/A_2
    U_2 = average_over_area(meanu_a, y_R, z_H)
#    U_2 = 1.0
    power_d = (U*ptot1 - U_2*ptot2)*(rho*A_2)
    eta2 = 0.5*rho*0.255/power_d
    print("eta_II =", eta2)
    print("A1/A2 =", A_1/A_2)
    print("p_2 =", p_2/rho)
    print("U_2 =", U_2)
    print("ptot_2 =", ptot2)
    print("Shaft power output:", 0.5*rho*1.0*0.255*1.0**3)
    print("Kinetic power dissipation:", power_d_k)
    print("Static power dissipation:", power_d_p)
    print("C_P/C_D =", 0.255/cd_meas)
   
def calc_re():
    a = np.load('Processed/a.npy')    
    U = 1.0
    tsr = np.arange(0.1, 3.2, 0.1)
    a = a[0:len(tsr)]
    U = (1-a)*U
    c = 0.14
    nu = 10**(-6)
    remed = U*tsr*c/nu
    remax = (U*tsr+U)*c/nu
    remin = np.abs((U*tsr-U)*c/nu)
    plt.close('all')
    plt.plot(tsr, remed/10**5, '-.k')
    plt.hold(True)
    plt.plot(tsr, remax/10**5, 'k', label=r'$Re_c$ bounds')
    plt.plot(tsr, remin/10**5, 'k')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$Re_c$ $(\times 10^5)$')
#    plt.legend(loc=4)
    styleplot()
    
def plot_torque_ripple():
    i = range(31)
    torque_ripple = np.load('Processed/torque_ripple.npy')
    tsr = np.load('Processed/tsr.npy')
    plt.plot(tsr[i], torque_ripple[i], '-ok', markerfacecolor = 'none')
    plt.xlabel(r'$\lambda$', labelpad=20)
    plt.ylabel(r'Torque ripple')
    styleplot()  
    plt.show()
    print("Torque ripple at TSR =", str(tsr[12])+":", torque_ripple[12])
    
def export_data():
    """Export processed data to text file."""
    import datetime
    datestring = datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")
    rev = 0
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
    metadata = ["Processed performance data from UNH RVAT experiment",
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
                "Bachant and Wosnik (2013) 'Performance and Wake Measurements for a Vertical Axis Turbine at Moderate Reynolds Number'", 
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
    with open('Processed/csv/unh-rvat-perf-2013-03-rev'+str(rev)+".txt",'wb') \
    as csvfile:
        fwriter = csv.writer(csvfile, delimiter='\t')
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
    plt.close("all")    
    plotsinglerun(100, perf=False, wake=False, autocorr=True)
#    plot_vel_spec(y_R=-0.1, z_H=0, tsr=1.9)
#    batchperf()
#    batchwake()
    sp = 'C:/Users/Pete/Google Drive/Research/Papers/JOT VAT near-wake/Figures/'
#    plotperf(save=True, savepath=sp)
#    plotwake(["fstrength"], save=False, savepath=sp)
    
if __name__ == "__main__":
    ts = time.time()
    main()
    te = time.time()
    print("Elapsed time:", te-ts, "s")