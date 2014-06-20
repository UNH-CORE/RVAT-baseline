# -*- coding: utf-8 -*-
"""
Created on Fri May 30 00:34:48 2014

@author: Pete
"""
from __future__ import division, print_function 
from processing import *

def setpltparams():
    font = {"family" : "serif", 
            "serif" : "cmr10",
            "sans-serif" : "cmr10",
            "size" : 23}
    lines = {"markersize" : 9, "markeredgewidth" : 0.9}
    legend = {"numpoints" : 1, "fontsize" : "small"}
    matplotlib.rc("text", usetex=True)
    matplotlib.rc("font", **font)
    matplotlib.rc("lines", **lines)
    matplotlib.rc("legend", **legend)
    matplotlib.rc("xtick", **{"major.pad":12})

def styleplot():
    plt.grid(True)
    plt.tight_layout()

def plotsinglerun(run, perf=True, wake=False, autocorr=False, save=False,
                  xaxis="time"):
    t1 = 13
    t2 = 30
    t2t = 30
    t, angle, Ttrans, Tarm, drag, rpm, tsr = loadtdms(run)
    t2, nrevs, nbp = find_t2(t, angle, t1, t2t)
    meantsr, std_tsr = calcstats(tsr, t1, t2, 2000)
    omega = meantsr*U/R
    blade_period = 2*np.pi/omega/3
    vecloaded = False
    if perf:
        cp = (Ttrans+0.5)*tsr/0.5/500.0
        ct = cp/tsr
        cd_ts = drag/500.0
        max_cd = np.max(cd_ts[t1*2000:t2*2000])
        min_cd = np.min(cd_ts[t1*2000:t2*2000])
        max_torque = np.max(Ttrans[t1*2000:t2*2000])
        min_torque = np.min(Ttrans[t1*2000:t2*2000])
        cd, std_cd = calcstats(cd_ts, t1, t2, 2000)
        meancp, std_cp = calcstats(cp, t1, t2, 2000)
        meanT, stdT = calcstats(Ttrans, t1, t2, 2000)
        meanrpm, std_rpm = calcstats(rpm, t1, t2, 2000)
        if xaxis == "time":
            x = t
            xlabel = r"$t$(s)"
        else:
            x = angle
            xlabel = r"Turbine rotation (deg)"
        plt.plot(x, cp, "k", label="Power")
        plt.plot(x, cd_ts, "r", label="Drag")
        plt.plot(x, ct, "--b", label="Torque")
        plt.xlabel(xlabel)
        plt.ylabel(r"Nondimensional coefficient")
        plt.legend(loc="lower left")
        if xaxis == "time":
            plot_vertical_lines([t1, t2t], color="black")
            plot_vertical_lines([t2], color="gray")
        styleplot()
        print("TSR =", meantsr, "; C_P =", meancp)
        print("Number of revolutions:", nrevs)
        print("(1/2)(max_cd-min_cd)/cd:", np.abs(max_cd-min_cd)*0.5/cd)
        print("(1/2)(Tmax-Tmin)/Tmean:", np.abs(max_torque-min_torque)*0.5/meanT)
    if wake:
        tv, u, v, w = loadvec(run)
        vecloaded = True
        angle = decimate(angle, 10)
        meanu, stdu = calcstats(u, t1, t2, 200)
        plt.figure()
        plt.plot(tv, u, "k")
        plt.xlabel("$t$ (s)")
        plt.ylabel("$u$ (m/s)")
        plt.vlines([t1, t2],np.min(u),np.max(u),
                   color="r",linestyles="dashed")
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
        # Find first zero crossing
#        i = np.where(np.diff(np.sign(rho)))[0][0]
#        int_time = np.trapz(rho[:i], x=tau[:i])
        int_time = np.trapz(rho, tau)
        print("Integral timescale =", int_time, "s")
        plt.figure()
        plt.plot(tau, rho)
        plot_vertical_lines([blade_period, blade_period*3], ymaxscale=1)
        plt.xlabel("Lag (s)")
        plt.ylabel("Autocorrelation coefficient")
        styleplot()
    plt.show()

def plotvelspec(y_R=0, z_H=0, tsr=1.9, newfig=True, show=False):
    """Plots the velocity spectrum for a single run."""
    # Find index for the desired parameters
    i = find_run_ind(y_R, z_H, tsr)
    print("Plotting spectra from run", i+1)
    t1 = 13
    t2 = pd.read_csv("Processed/processed.csv")["t2"][i]
    t, u, v, w = loadvec(i+1) # Run name is index + 1
    v_seg = v[200*t1:200*t2] - np.mean(v[200*t1:200*t2])
    f, spec = psd(t, v_seg, window="Hanning")
    f_turbine = tsr*U/R/(2*np.pi)
    f_blade = f_turbine*3
    # Find maximum frequency and its relative strength
    f_max = f[np.where(spec==np.max(spec))[0][0]]
    strength = np.max(spec)/np.var(v_seg)*(f[1] - f[0])
    print("Strongest frequency f/f_turbine:", f_max/f_turbine)
    print("Relative strength:", strength)
    # Calculate shaft shedding frequency
    St = 0.19 # Approximate for Re_d = 1e5
    f_cyl = St*U/d_shaft
    if newfig:
        plt.figure()
    plt.loglog(f/f_turbine, spec, "k")
    plt.xlim((0, 50))
    plt.xlabel(r"$f/f_{\mathrm{turbine}}$")
    plt.ylabel(r"Spectral density")
    # Should the spectrum be normalized?
    f_line = np.linspace(10,40)
    spec_line = f_line**(-5./3)*0.1
    plt.hold(True)
    plt.loglog(f_line, spec_line)
    plt.ylim((10**-9, 1))
    plot_vertical_lines([1, 3, 6, 9])
    styleplot()
    plt.grid()
    if show:
        plt.show()
        
def plotperfspec(y_R=0, z_H=0, tsr=1.9, newfig=True, show=False):
    """Plots the performance spectra for a single run."""
    # Find index for the desired parameters
    i = find_run_ind(y_R, z_H, tsr)
    print("Plotting spectra from run", i+1)
    t1 = 13
    t2 = pd.read_csv("Processed/processed.csv")["t2"][i]
    t, angle, Ttrans, Tarm, drag, rpm, tsr_ts = loadtdms(i+1) # Run name is index + 1
    torque = Tarm/(0.5*rho*A_t*R*U**2)
    torque_seg = torque[2000*t1:2000*t2] - np.mean(torque[2000*t1:2000*t2])
    f, spec = psd(t, torque_seg, window="Hanning")
    f_turbine = tsr*U/R/(2*np.pi)
    f_blade = f_turbine*3
    # Find maximum frequency and its relative strength
    f_max = f[np.where(spec==np.max(spec))[0][0]]
    strength = np.max(spec)/np.var(torque_seg)*(f[1] - f[0])
    print("Strongest frequency f/f_turbine:", f_max/f_turbine)
    print("Spectral concentration:", strength)
    if newfig:
        plt.figure()
    plt.loglog(f/f_turbine, spec, "k")
    plt.xlim((0, 50))
    plt.xlabel(r"$f/f_{\mathrm{turbine}}$")
    plt.ylabel(r"Spectral density")
    # Should the spectrum be normalized?
    plot_vertical_lines([1, 3, 6, 9])
    styleplot()
    plt.grid()
    if show:
        plt.show()
        
def plotmultispec(save=False, savepath="", savetype=".pdf"):
    """Creates a 1x3 plot for spectra of torque coefficient and cross-stream
    velocity spectra at two locations."""
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 3, 1)
    plotperfspec(y_R=-1, z_H=0.25, tsr=1.9, newfig=False)
    plt.title("(a)", fontsize=20)
    plt.subplot(1, 3, 2)
    plotvelspec(y_R=-1, z_H=0.25, tsr=1.9, newfig=False)
    plt.title("(b)", fontsize=20)
    plt.ylabel("")
    plt.subplot(1, 3, 3)
    plotvelspec(y_R=1.5, z_H=0.25, tsr=1.9, newfig=False)
    plt.title("(c)", fontsize=20)
    plt.ylabel("")
    plt.annotate(r"$f^{-5/3}$", xy=(12, 0.5e-2), fontsize=16)
    plt.tight_layout()
    if save:
        plt.savefig(savepath + "/multispec" + savetype)
    plt.show()
    
def plot_vertical_lines(xlist, ymaxscale=1, color="gray"):
    if not isinstance(xlist, list):
        x = [x]
    ymin = plt.axis()[2]
    ymax = plt.axis()[3]*ymaxscale
    for x in xlist:
        plt.vlines(x, ymin, ymax,
                   color=color, linestyles="dashed")
    plt.ylim((ymin, ymax))
    
def plotvelhist(run):
    """Plots the velocity histogram for a given run."""
    i = run - 1 # Run indexing starts from 1!
    t1 = 13
    t2 = pd.read_csv("Processed/processed.csv")["t2"][i]
    t, u, v, w = loadvec(run)
    u = u[t1*200:t2*200]
    plt.figure()
    plt.hist(u-u.mean(), bins=50, histtype="step", color="k", normed=True)
    plt.xlabel(r"$u-U$")
    plt.ylabel("Samples (normalized)")
    styleplot()
    plt.grid(False)
    plt.show()
    
def plotwake(plotlist, save=False, savepath=None, savetype=".pdf",
             print_analysis=False):
    if not isinstance(plotlist, list):
        plotlist = [plotlist]
    z_H = np.arange(0, 0.75, 0.125)
    y_R = np.hstack([-3.,-2.75,-2.5,-2.25,-2.,-1.8,np.arange(-1.6,0.1,0.1)])
    y_R = np.hstack([y_R, -np.flipud(y_R[0:-1])])
    y_R = np.round(y_R, decimals=4)
    # Load processed data
    df = pd.read_csv("Processed/processed.csv")
    df["k"] = 0.5*(df.stdu**2 + df.stdv**2 + df.stdw**2)
    df["meank"] = 0.5*(df.meanu**2 + df.meanv**2 + df.meanw**2)
    df["kbar"] = df.meank + df.k
    # Create empty 2D arrays for contour plots, etc.
    grdata = {}
    grdims = (len(z_H), len(y_R))
    quantities = ["meanu", "meanv", "meanw", "stdu", "stdv", "stdw",
                  "meanupvp", "meanupwp", "meanvpwp", "meanupup", "meanvpvp", 
                  "meanwpwp", "meanuu", "vectemp", "fpeak_u", "fstrength_u", 
                  "fpeak_v", "fstrength_v", "fpeak_w", "fstrength_w", "k",
                  "meank", "kbar"]
    for q in quantities:
        grdata[q] = np.zeros(grdims)
    # Populate 2D arrays for velocity fields
    for n in range(len(z_H)):
        runs = getruns(z_H[n], tsr=1.9)
        i = [run - 1 for run in runs]
        for q in quantities:
            grdata[q][n,:] = df[q][i]
    def turb_lines():
        plt.hlines(0.5, -1, 1, linestyles="solid", linewidth=2)
        plt.vlines(-1, 0, 0.5, linestyles="solid", linewidth=2)
        plt.vlines(1, 0, 0.5, linestyles="solid", linewidth=2)
    def calc_meankturbtrans():
        z = H*z_H
        y = R*y_R
        ddy_uvU = np.zeros(grdims)
        ddz_uwU = np.zeros(grdims)
        ddy_vvV = np.zeros(grdims)
        ddz_vwV = np.zeros(grdims)
        ddy_vwW = np.zeros(grdims)
        ddz_wwW = np.zeros(grdims)
        meanu, meanv, meanw = grdata["meanu"], grdata["meanv"], grdata["meanw"]
        uv, vv, vw = grdata["meanupvp"], grdata["meanvpvp"], grdata["meanvpwp"]
        uw, ww = grdata["meanupwp"], grdata["meanwpwp"]
        for n in range(len(z)):
            ddy_uvU[n,:] = fdiff.second_order_diff((uv*meanu)[n,:], y)
            ddy_vvV[n,:] = fdiff.second_order_diff((vv*meanv)[n,:], y)
            ddy_vwW[n,:] = fdiff.second_order_diff((vw*meanw)[n,:], y)
        for n in range(len(y)):
            ddz_uwU[:,n] = fdiff.second_order_diff((uw*meanu)[:,n], z)
            ddz_vwV[:,n] = fdiff.second_order_diff((vw*meanv)[:,n], z)
            ddz_wwW[:,n] = fdiff.second_order_diff((ww*meanw)[:,n], z)
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
        for n in range(len(z)):
            dUdy[n,:] = fdiff.second_order_diff(meanu_a[n,:], y)
            dVdy[n,:] = fdiff.second_order_diff(meanv_a[n,:], y)
            dWdy[n,:] = fdiff.second_order_diff(meanw_a[n,:], y)
        for n in range(len(y)):
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
        for n in range(len(z)):
            dKdy[n,:] = fdiff.second_order_diff(meank_a[n,:], y)
        for n in range(len(y)):
            dKdz[:,n] = fdiff.second_order_diff(meank_a[:,n], z)
        return dKdy, dKdz
    if "meanucont" in plotlist or "all" in plotlist:
        # Plot contours of mean streamwise velocity
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, meanu_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.3)
        cb.set_label(r"$\overline{u}/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/meanucont"+savetype)
    if "v-wquiver" in plotlist or "all" in plotlist:
        # Make quiver plot of v and w velocities
        plt.figure(figsize=(10,5))
        Q = plt.quiver(y_R, z_H,meanv_a, meanw_a, angles="xy")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.2, 0.1, r"$0.1$ m/s",
                   labelpos="E",
                   coordinates="figure",
                   fontproperties={"size": "small"})
        plt.tight_layout()
        plt.hlines(0.5, -1, 1, linestyles="solid", colors="r",
                   linewidth=2)
        plt.vlines(-1, -0.2, 0.5, linestyles="solid", colors="r",
                   linewidth=2)
        plt.vlines(1, -0.2, 0.5, linestyles="solid", colors="r",
                   linewidth=2)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/v-wquiver"+savetype)
    if "meanu_2tsrs" in plotlist:
        # Plot mean velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanu[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanu[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{u}/U_\infty$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/meanu_2tsrs"+savetype)
    if "stdu_2tsrs" in plotlist:
        # Plot stdu velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdu[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdu[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\sigma_u/U_\infty$")
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+"/stdu_2tsrs"+savetype)
    if "uw_2tsrs" in plotlist:
        # Plot uw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuw[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuw[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{u'w'}$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/uw_2tsrs"+savetype)
    if "meanuvstsr" in plotlist:
        # Plot mean velocity components vs TSR
        tsr = np.load("Processed/tsr.npy")
        runs = range(1,32)
        ind = [run-1 for run in runs]
        plt.figure()
        plt.plot(tsr[ind], meanu[ind], "-ok", markerfacecolor="none")
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$\overline{u}/U_\infty$")
        plt.hold(True)
        plt.plot(tsr[ind], meanv[ind], "-sk", markerfacecolor="none")
        plt.plot(tsr[ind], meanw[ind], "-^k", markerfacecolor="none")
        runs = range(347,378)
        ind = [run-1 for run in runs]
        plt.plot(tsr[ind], meanu[ind], "-ok", markerfacecolor="k",
             label=r"$\overline{u}$")
        plt.plot(tsr[ind], meanv[ind], "-sk", markerfacecolor="k",
             label=r"$\overline{v}$")
        plt.plot(tsr[ind], meanw[ind], "-^k", markerfacecolor="k",
             label=r"$\overline{w}$")
        plt.legend(ncol=3)
        styleplot()
        if save:
            plt.savefig(savepath+"/meanuvstsr"+savetype)
    if "stducont" in plotlist:
        # Plot contours of streamwise turbulence intensity
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, stdu_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.3)
        cb2.set_label(r"$\sigma_u/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/stducont"+savetype)
    if "uvcont" in plotlist:
        # Plot contours of uv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanuv_a, 15, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$\overline{u'v'}/U_\infty^2$")
        cb2.set_ticks(np.arange(-0.02, 0.025, 0.005), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/uvcont"+savetype)
    if "meanw_2tsrs" in plotlist:
        # Plot mean vertical velocity profiles at two TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanw[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanw[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{w}/U_\infty$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/meanw_2tsrs"+savetype)
    if "meanv_2tsrs" in plotlist:
        # Plot mean cross stream velocity profiles at two TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanv[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanv[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{v}/U_\infty$")
        plt.ylim(-0.1501, 0.1001)
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/meanv_2tsrs"+savetype)
    if "stdv_2tsrs" in plotlist:
        # Plot stdv velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdv[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdv[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\sigma_v/U_\infty$")
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+"/stdv_2tsrs"+savetype)
    if "stdw_2tsrs" in plotlist:
        # Plot stdw velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdw[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, stdw[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\sigma_w/U_\infty$")
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+"/stdw_2tsrs"+savetype)
    if "uv_2tsrs" in plotlist:
        # Plot uv Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuv[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanuv[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{u'v'}$")
        plt.legend()
        styleplot()
        if save:
            plt.savefig(savepath+"/uv_2tsrs"+savetype)
    if "kcont" in plotlist or "all" in plotlist:
        # Plot contours of k
        plt.figure(figsize=(10,5))
        csphi = plt.contourf(y_R, z_H, k_a/(0.5*1.0**2), 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cbphi = plt.colorbar(csphi, shrink=1, extend="both", 
                             orientation="horizontal", pad=0.26)
        cbphi.set_label(r"$k/\frac{1}{2}U_\infty^2$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/kcont"+savetype)
    if "meankcont" in plotlist:
        # Plot contours of k
        plt.figure(figsize=(10,5))
        csphi = plt.contourf(y_R, z_H, meank_a/(0.5*1**2), 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cbphi = plt.colorbar(csphi, shrink=1, extend="both", 
                             orientation="horizontal", pad=0.26)
        cbphi.set_label(r"$K/\frac{1}{2}U_\infty^2$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/meankcont"+savetype)
    if "meanvcont" in plotlist:
        # Plot contours of meanv
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, meanv_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cbmv.set_label(r"$\overline{v}/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/meanvcont"+savetype)
    if "stdvcont" in plotlist:
        # Plot contours of stdv
        plt.figure(figsize=(10,5))
        cstdv = plt.contourf(y_R, z_H, stdv_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbstdv = plt.colorbar(cstdv, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cbstdv.set_label(r"$\sigma_v/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/stdvcont"+savetype)
    if "meanwcont" in plotlist:
        # Plot contours of meanw
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, meanw_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cbmv.set_label(r"$\overline{w}/U_{\infty}$")
        turb_lines()
        if save:
            plt.savefig(savepath+"/meanwcont"+savetype)
    if "stdwcont" in plotlist:
        # Plot contours of stdw
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, stdw_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cbmv.set_label(r"$\sigma_w/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/stdwcont"+savetype)
    if "vw_2tsrs" in plotlist:
        # Plot vw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvw[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvw[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{v'w'}$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/vw_2tsrs"+savetype)
    if "vwcont" in plotlist:
        # Plot contours of vw Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanvw_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.3)
        cb2.set_label(r"$\overline{v'w'}/U_\infty^2$")
        cb2.set_ticks(np.linspace(-.008,.006,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/vwcont"+savetype)
    if "uwcont" in plotlist:
        # Plot contours of vw Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanuw_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$\overline{u'w'}/U_\infty^2$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/uwcont"+savetype)
    if "vvcont" in plotlist:
        # Plot contours of vv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanvv_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cb2.set_label(r"$\overline{v'v'}$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/vvcont"+savetype)
    if "wwcont" in plotlist:
        # Plot contours of vv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanww_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cb2.set_label(r"$\overline{w'w'}$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/wwcont"+savetype)
    if "uucont" in plotlist:
        # Plot contours of uu Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, meanuu_a, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.2)
        cb2.set_label(r"$\overline{u"u"}$")
        cb2.set_ticks(np.linspace(0,.108,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/uucont"+savetype)
    if "vv_2tsrs" in plotlist:
        # Plot vw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvv[ind], "-ok", markerfacecolor="none", 
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, meanvv[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{v'v'}$")
        plt.legend()
        styleplot()
        if save:
            plt.savefig(savepath+"/vv_2tsrs"+savetype)
    if "fpeak_u" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fpeak_u_a, cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$f_{\mathrm{peak}}/f_{\mathrm{turbine}}$")
        cb2.set_ticks(np.linspace(0,10,11), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/fpeak_u"+savetype)
    if "fstrength_u" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fstrength_u_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$S_{\max}/\sigma^2_u$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/fstrength_u"+savetype)
    if "fpeak_v" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fpeak_v_a, cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$f_{\mathrm{peak}}/f_{\mathrm{turbine}}$")
        cb2.set_ticks(np.linspace(0,10,11), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/fpeak_v"+savetype)
    if "fstrength_v" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fstrength_v_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$\Psi$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/fstrength_v"+savetype)
    if "fpeak_w" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fpeak_w_a, cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$f_{\mathrm{peak}}/f_{\mathrm{turbine}}$")
        cb2.set_ticks(np.linspace(0,10,11), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/fpeak_w"+savetype)
    if "fstrength_w" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, fstrength_w_a, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both", 
                           orientation="horizontal", pad=0.26)
        cb2.set_label(r"$S_{\max}/\sigma^2_w$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"fstrength_w"+savetype)
    # Plot estimate for production of turbulence kinetic energy
    if "kprod" in plotlist:
        calc_meanvelgrad()
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, kprod, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.26)
        cb.set_label(r"$-\overline{u_i' u_j'}\frac{\partial U_i}{\partial x_j}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/kprod"+savetype)
    if "meankadv" in plotlist:
        calc_meankgrad()
        # Make quiver plot of K advection
        plt.figure(figsize=(10,5))
        plt.hlines(0.5, -1, 1, linestyles="solid", colors="r",
                   linewidth=3)
        plt.vlines(-1, -0.2, 0.5, linestyles="solid", colors="r",
                   linewidth=3)
        plt.vlines(1, -0.2, 0.5, linestyles="solid", colors="r",
                   linewidth=3)
        Q = plt.quiver(y_R, z_H, meanv_a/meanu_a*dKdy, meanw_a/meanu_a*dKdz, 
                       scale=4, angles="xy")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.1, 0.2, r"$0.2 \mathrm{\, m/s^2}$",
                      labelpos="E",
                      coordinates="figure",
                      fontproperties={"size": "small"})
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.grid(False)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/meankadv"+savetype)
    if "Kturbtrans" in plotlist:
        calc_meankturbtrans()
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, tty, 20, cmap=plt.cm.coolwarm,
                          levels=np.linspace(-0.08, 0.08, 21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.26)
        cb.set_label(r"$-\frac{1}{2}\frac{\partial}{\partial x_j}\overline{u_i' u_j'} U_i$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/Kturbtrans"+savetype)
    if "meancomboquiv" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,6))
        # Add contours of mean velocity
        cs = plt.contourf(y_R, z_H, meanu_a, 20, cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cb.set_label(r"$U/U_{\infty}$")
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(y_R, z_H, meanv_a, meanw_a, angles="xy", width=0.0022)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.3, 0.1, r"$0.1 U_\infty$",
                   labelpos="E",
                   coordinates="figure",
                   fontproperties={"size": "small"})
        plt.hlines(0.5, -1, 1, linestyles="solid", colors="gray",
                   linewidth=3)
        plt.vlines(-1, -0.2, 0.5, linestyles="solid", colors="gray",
                   linewidth=3)
        plt.vlines(1, -0.2, 0.5, linestyles="solid", colors="gray",
                   linewidth=3)
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/meancomboquiv"+savetype)
    if "xvorticity" in plotlist:
        z = 1.0*z_H
        y = R*y_R
        dWdy = np.zeros(np.shape(meanu_a))
        dVdz = np.zeros(np.shape(meanu_a))
        for n in range(len(z)):
            dWdy[n,:] = fdiff.second_order_diff(meanw_a[n,:], y)
        for n in range(len(y)):
            dVdz[:,n] = fdiff.second_order_diff(meanv_a[:,n], z)
        # Make quiver plot of K advection
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, dWdy-dVdz, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.26)
        cb.set_label(r"$\Omega_x$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/xvorticity"+savetype)
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
        quantities = [average_over_area(-2*meanv_a/meanu_a*dKdy/(0.5*U**2)/D, y_R, z_H), 
                      average_over_area(-2*meanw_a/meanu_a*dKdz/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*tty/meanu_a/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*ttz/meanu_a/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*kprod/meanu_a/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*meandiss/meanu_a/(0.5*U**2)/D, y_R, z_H)]
        ax = plt.gca()
        ax.bar(range(len(names)), quantities, color="k", width=0.5)
        ax.set_xticks(np.arange(len(names))+0.25)
        ax.set_xticklabels(names)
        plt.hlines(0, 0, len(names), color="gray")
        plt.ylabel(r"$\frac{K \, \mathrm{ transport}}{UDK_\infty}$")
#        ax.annotate(r"$\mathrm{Total} = " \
#                    + str(np.round(np.sum(quantities), decimals=4)) + "$", 
#                    xy=(0, 0), xytext=(0.75, 0.82), 
#                    xycoords="figure fraction", fontsize=18)
        styleplot()
        plt.grid(False)
        if print_analysis:
            print("K recovery rate (%/D) =", 
                  2*np.sum(quantities)/(0.5*U**2)/D*100)
        if save:
            plt.savefig(savepath+"/Kbargraph"+savetype)
    plt.show(block=False)
    if print_analysis:
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
    
def plot_torque_ripple():
    i = range(31)
    torque_ripple = np.load("Processed/torque_ripple.npy")
    tsr = np.load("Processed/tsr.npy")
    plt.plot(tsr[i], torque_ripple[i], "-ok", markerfacecolor = "none")
    plt.xlabel(r"$\lambda$", labelpad=20)
    plt.ylabel(r"Torque ripple")
    styleplot()  
    plt.show()
    print("Torque ripple at TSR =", str(tsr[12])+":", torque_ripple[12])
    
def plot_Re_c():
    a = np.load("Processed/a.npy")    
    U = 1.0
    tsr = np.arange(0.1, 3.2, 0.1)
    a = a[0:len(tsr)]
    U = (1-a)*U
    c = 0.14
    nu = 10**(-6)
    remed = U*tsr*c/nu
    remax = (U*tsr+U)*c/nu
    remin = np.abs((U*tsr-U)*c/nu)
    plt.close("all")
    plt.plot(tsr, remed/10**5, "-.k")
    plt.hold(True)
    plt.plot(tsr, remax/10**5, "k", label=r"$Re_c$ bounds")
    plt.plot(tsr, remin/10**5, "k")
    plt.xlabel(r"$\lambda$")
    plt.ylabel(r"$Re_c$ $(\times 10^5)$")
#    plt.legend(loc=4)
    styleplot()

def plotperf(plotlist=["cp", "cd"],
             subplots=True, save=False, savepath="", savetype=".pdf"):
    i = np.arange(31)
    data = pd.read_csv("Processed/processed.csv")
    cp = data.cp
    cd = data.cd
    tsr = data.tsr
    eta2 = data.eta2
    if not subplots:
        # Exergy efficiency plot
        plt.figure()
        plt.plot(tsr[i], eta2[i], "-ok", markerfacecolor = "none")
        plt.xlabel(r"$\lambda$", labelpad=20)
        plt.ylabel(r"$\eta_{II}$")
        styleplot()
        if save:
            plt.savefig(savepath+"/eta2"+savetype)
        # Power coefficient plot
        plt.figure()
        plt.plot(tsr[i], cp[i], "-ok", markerfacecolor = "none")
        plt.xlabel(r"$\lambda$", labelpad=20)
        plt.ylabel(r"$C_P$")
        styleplot()
        if save:
            plt.savefig(savepath+"/cpvstsr"+savetype)
        # Drag coefficient plot
        plt.figure()
        plt.plot(tsr[i], cd[i], "-ok", markerfacecolor = "none")
        plt.xlabel(r"$\lambda$", labelpad=20)
        plt.ylabel(r"$C_D$")
        plt.ylim(0, 1.2001)
        styleplot()
        if save:
            plt.savefig(savepath+"/cdvstsr"+savetype)
        # Torque coefficient plot
        ct = cp/tsr
        plt.figure()
        plt.plot(tsr[i], ct[i], "-ok", markerfacecolor = "none")
        plt.xlabel(r"$\lambda$", labelpad=20)
        plt.ylabel(r"$C_T$")
        styleplot()
        if save:
            plt.savefig(savepath+"/ctvstsr"+savetype)
    else:
        plt.figure(figsize=(12,5))
        plt.subplot(121)
        plt.plot(tsr[i], cp[i], "-ok", markerfacecolor = "none")
        plt.xlabel(r"$\lambda$", labelpad=20)
        plt.ylabel(r"$C_P$")
        plt.grid()
        plt.subplot(122)
        plt.plot(tsr[i], cd[i], "-ok", markerfacecolor = "none")
        plt.xlabel(r"$\lambda$", labelpad=20)
        plt.ylabel(r"$C_D$")
        plt.ylim(0, 1.2001)
        styleplot()
        if save:
            plt.savefig(savepath+"/perf"+savetype)
    print("At tsr = 1.9, C_P =", cp[np.where(np.round(tsr, decimals=2)==1.9)[0][0]],
          "; C_D =", cd[np.where(np.round(tsr, decimals=2)==1.9)[0][0]])
    plt.show()
        
def plotperf_periodic():
    i = range(31)
    d = pd.read_csv("Processed/processed.csv")
    plt.figure()
    plt.plot(d.tsr[i], d.amp_cd[i])
    styleplot()
    plt.figure()
    plt.plot(d.tsr[i], d.phase_cd[i])
        
def main():
    setpltparams()
    plt.close("all")
    p = "Google Drive/Research/Papers/JOT CFT near-wake/Figures"
    if "linux" in sys.platform:
        p = "/home/pete/" + p
    elif "win" in sys.platform:
        p = "C:/Users/Pete/" + p
        
#    plotsinglerun(111, perf=True, wake=False, autocorr=False, xaxis="angle")
#    plotvelspec(y_R=1.5, z_H=0.25, tsr=1.9, show=True)
#    plotperfspec(y_R=1.5, z_H=0.25, tsr=1.9, show=True)
#    plotperf(subplots=True, save=False, savepath=p)
    plotwake(["meanu"], save=False, savepath=p,
             print_analysis=True)
#    plotmultispec(save=False, savepath=p)
#    plotperf_periodic()
#    plotvelhist(5)
        
if __name__ == "__main__":
    main()