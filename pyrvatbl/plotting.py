# -*- coding: utf-8 -*-
"""This module contains plotting functions."""

from __future__ import division, print_function
from .processing import *
import scipy.stats
from pxl.styleplot import set_sns


def setpltparams(seaborn=True, fontsize=18, latex=True):
    if seaborn:
        set_sns()
    else:
        if latex:
            font = {"family" : "serif", "serif" : "cmr10", "size" : fontsize}
        else:
            font = {"size" : fontsize}
        lines = {"markersize" : 9, "markeredgewidth" : 1,
                 "linewidth" : 1.2}
        legend = {"numpoints" : 1, "fontsize" : "small"}
        matplotlib.rc("text", usetex=latex)
        matplotlib.rc("font", **font)
        matplotlib.rc("lines", **lines)
        matplotlib.rc("legend", **legend)
        matplotlib.rc("xtick", **{"major.pad":6})


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
        print("TSR =", meantsr, "\nC_P =", meancp)
        print("Number of revolutions: {:d}".format(int(nrevs)))
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
        tau, rho = timeseries.autocorr_coeff(u, t, 0, 6.0)
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


def plotvelspec(y_R=0, z_H=0, tsr=1.9, newfig=True, n_band_average=1,
                plot_conf_int=False, verbose=False):
    """Plots the velocity spectrum for a single run."""
    # Find index for the desired parameters
    i = find_run_ind(y_R, z_H, tsr)
    if verbose:
        print("Plotting spectra from run {}".format(i+1))
    t1 = 13
    t2 = pd.read_csv("Data/Processed/processed.csv")["t2"][i]
    t, u, v, w = loadvec(i+1) # Run name is index + 1
    v_seg = v[200*t1:200*t2] - np.mean(v[200*t1:200*t2])
    f, spec = psd(t, v_seg, window="Hanning", n_band_average=n_band_average)
    f_turbine = tsr*U/R/(2*np.pi)
    # Find maximum frequency and its relative strength
    f_max = f[np.where(spec==np.max(spec))[0][0]]
    strength = np.max(spec)/np.var(v_seg)*(f[1] - f[0])
    if verbose:
        print("Strongest frequency f/f_turbine: {:.3f}".format(f_max/f_turbine))
        print("Spectral concentration: {:.3f}".format(strength))
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
    spec_line = f_line**(-5./3)*0.5
    plt.hold(True)
    plt.loglog(f_line, spec_line, "gray")
    plt.ylim((10**-9, 1))
    plot_vertical_lines([1, 3, 6, 9])
    if plot_conf_int:
        dof = n_band_average*2
        chi2 = scipy.stats.chi2.interval(alpha=0.95, df=dof)
        y1 = dof*spec/chi2[1]
        y2 = dof*spec/chi2[0]
        plt.fill_between(f/f_turbine, y1, y2, facecolor="lightgray", alpha=0.2)
    plt.tight_layout()


def plotperfspec(y_R=0, z_H=0, tsr=1.9, newfig=True, verbose=False,
                 n_band_average=1, plot_conf_int=False):
    """Plots the performance spectra for a single run."""
    # Find index for the desired parameters
    i = find_run_ind(y_R, z_H, tsr)
    if verbose:
        print("Plotting spectra from run {}".format(i+1))
    t1 = 13
    t2 = pd.read_csv("Data/Processed/processed.csv")["t2"][i]
    t, angle, Ttrans, Tarm, drag, rpm, tsr_ts = loadtdms(i + 1)
    torque = Tarm/(0.5*rho*A_t*R*U**2)
    torque_seg = torque[2000*t1:2000*t2] - np.mean(torque[2000*t1:2000*t2])
    f, spec = psd(t, torque_seg, window="Hanning",
                  n_band_average=n_band_average)
    f_turbine = tsr*U/R/(2*np.pi)
    # Find maximum frequency and its relative strength
    f_max = f[np.where(spec==np.max(spec))[0][0]]
    strength = np.max(spec)/np.var(torque_seg)*(f[1] - f[0])
    if verbose:
        print("Strongest frequency f/f_turbine: {:.3f}".format(f_max/f_turbine))
        print("Spectral concentration: {:.3f}".format(strength))
    if newfig:
        plt.figure()
    plt.loglog(f/f_turbine, spec, "k")
    plt.xlim((0, 50))
    plt.xlabel(r"$f/f_{\mathrm{turbine}}$")
    plt.ylabel(r"Spectral density")
    # Should the spectrum be normalized?
    if plot_conf_int:
        dof = n_band_average*2
        chi2 = scipy.stats.chi2.interval(alpha=0.95, df=dof)
        y1 = dof*spec/chi2[1]
        y2 = dof*spec/chi2[0]
        plt.fill_between(f/f_turbine, y1, y2, facecolor="lightgray", alpha=0.2)
    plot_vertical_lines([1, 3, 6, 9])
    plt.tight_layout()


def plotmultispec(save=False, savepath="Figures", savetype=".pdf",
                  n_band_average=5, plot_conf_int=False, verbose=False):
    """Creates a 1x3 plot for spectra of torque coefficient and cross-stream
    velocity spectra at two locations."""
    plt.figure(figsize=(7.5, 3.75/10*7.5))
    plt.subplot(1, 3, 1)
    plotperfspec(y_R=-1, z_H=0.25, tsr=1.9, newfig=False, verbose=verbose,
                 n_band_average=n_band_average, plot_conf_int=plot_conf_int)
    plt.title("(a)")
    plt.subplot(1, 3, 2)
    plotvelspec(y_R=-1, z_H=0.25, tsr=1.9, newfig=False, verbose=verbose,
                n_band_average=n_band_average, plot_conf_int=plot_conf_int)
    plt.title("(b)")
    plt.ylabel("")
    plt.annotate(r"$f^{-5/3}$", xy=(12, 1.5e-2), fontsize=12)
    plt.subplot(1, 3, 3)
    plotvelspec(y_R=1.5, z_H=0.25, tsr=1.9, newfig=False, verbose=verbose,
                n_band_average=n_band_average, plot_conf_int=plot_conf_int)
    plt.title("(c)")
    plt.ylabel("")
    plt.annotate(r"$f^{-5/3}$", xy=(12, 1.5e-2), fontsize=12)
    plt.tight_layout()
    if save:
        plt.savefig(savepath + "/multispec" + savetype)


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
    t2 = pd.read_csv("Data/Processed/processed.csv")["t2"][i]
    t, u, v, w = loadvec(run)
    u = u[t1*200:t2*200]
    plt.figure()
    plt.hist(u-u.mean(), bins=50, histtype="step", color="k", normed=True)
    plt.xlabel(r"$u-U$")
    plt.ylabel("Samples (normalized)")
    styleplot()
    plt.grid(False)


def plotwake(plotlist, scale=1, save=False, savepath="Figures",
             savetype=".pdf", print_analysis=False, barcolor="gray"):
    if not isinstance(plotlist, list):
        plotlist = [plotlist]
    figsize_horiz_contour = np.array((7.5, 4.5))*scale
    figsize_vertical_contour = np.array((7.5, 1.81))*scale
    figsize_vertical_quiver = np.array((7.5, 2.5))*scale
    horiz_cbar_pad = 0.17
    vertical_cbar_pad = 0.02
    vertical_cbar_shrink = 0.9
    # Load processed data
    df = pd.read_csv("Data/Processed/processed.csv")
    df["k"] = 0.5*(df.stdu**2 + df.stdv**2 + df.stdw**2)
    df["meank"] = 0.5*(df.meanu**2 + df.meanv**2 + df.meanw**2)
    df["kbar"] = df.meank + df.k
    # Create empty 2D arrays for contour plots, etc.
    quantities = ["meanu", "meanv", "meanw", "stdu", "stdv", "stdw",
                  "meanupvp", "meanupwp", "meanvpwp", "meanupup", "meanvpvp",
                  "meanwpwp", "meanuu", "vectemp", "fpeak_u", "fstrength_u",
                  "fpeak_v", "fstrength_v", "fpeak_w", "fstrength_w", "k",
                  "meank", "kbar"]
    # Create DataFrame pivoted for velocity field contour and quiver plots
    i = np.arange(31, 301)
    grdata = df[quantities + ["y/R", "z/H"]].iloc[i]
    grdata = grdata.pivot(index="z/H", columns="y/R")
    y_R = grdata.meanu.columns.values
    z_H = grdata.index.values
    z = H*z_H
    y = R*y_R
    grdims = (len(z_H), len(y_R))
    # Create some global variables from the grid data for cleaner syntax
    meanu, meanv, meanw = grdata["meanu"], grdata["meanv"], grdata["meanw"]
    uv, vv, vw = grdata["meanupvp"], grdata["meanvpvp"], grdata["meanvpwp"]
    uw, ww = grdata["meanupwp"], grdata["meanwpwp"]
    def turb_lines():
        plt.hlines(0.5, -1, 1, linestyles="solid", linewidth=2)
        plt.vlines(-1, 0, 0.5, linestyles="solid", linewidth=2)
        plt.vlines(1, 0, 0.5, linestyles="solid", linewidth=2)
    def calc_meankturbtrans():
        ddy_uvU = np.zeros(grdims)
        ddz_uwU = np.zeros(grdims)
        ddy_vvV = np.zeros(grdims)
        ddz_vwV = np.zeros(grdims)
        ddy_vwW = np.zeros(grdims)
        ddz_wwW = np.zeros(grdims)
        for n in range(len(z)):
            ddy_uvU[n,:] = fdiff.second_order_diff((uv*meanu).iloc[n,:], y)
            ddy_vvV[n,:] = fdiff.second_order_diff((vv*meanv).iloc[n,:], y)
            ddy_vwW[n,:] = fdiff.second_order_diff((vw*meanw).iloc[n,:], y)
        for n in range(len(y)):
            ddz_uwU[:,n] = fdiff.second_order_diff((uw*meanu).iloc[:,n], z)
            ddz_vwV[:,n] = fdiff.second_order_diff((vw*meanv).iloc[:,n], z)
            ddz_wwW[:,n] = fdiff.second_order_diff((ww*meanw).iloc[:,n], z)
        tt = -0.5*(ddy_uvU + ddz_uwU + ddy_vvV + ddz_vwV + ddy_vwW + ddz_wwW)
        tty = -0.5*(ddy_uvU + ddy_vvV + ddy_vwW) # Only ddy terms
        ttz = -0.5*(ddz_uwU + ddz_vwV + ddz_wwW) # Only ddz terms
        return tt, tty, ttz
    def calc_kprod_meandiss():
        dUdy = np.zeros(grdims)
        dUdz = np.zeros(grdims)
        dVdy = np.zeros(grdims)
        dVdz = np.zeros(grdims)
        dWdy = np.zeros(grdims)
        dWdz = np.zeros(grdims)
        for n in range(len(z)):
            dUdy[n,:] = fdiff.second_order_diff(meanu.iloc[n,:], y)
            dVdy[n,:] = fdiff.second_order_diff(meanv.iloc[n,:], y)
            dWdy[n,:] = fdiff.second_order_diff(meanw.iloc[n,:], y)
        for n in range(len(y)):
            dUdz[:,n] = fdiff.second_order_diff(meanu.iloc[:,n], z)
            dVdz[:,n] = fdiff.second_order_diff(meanv.iloc[:,n], z)
            dWdz[:,n] = fdiff.second_order_diff(meanw.iloc[:,n], z)
        kprod = uv*dUdy + uw*dUdz + vw*dVdz + vw*dWdy\
                + vv*dVdy + ww*dWdz
        meandiss = -2.0*nu*(dUdy**2 + dUdz**2 + dVdy**2 + dVdz**2 + dWdy**2 + dWdz**2)
        return kprod, meandiss
    def calc_meankgrad():
        z = H*z_H
        y = R*y_R
        dKdy = np.zeros(np.shape(meanu))
        dKdz = np.zeros(np.shape(meanu))
        for n in range(len(z)):
            dKdy[n,:] = fdiff.second_order_diff(grdata.meank.iloc[n,:], y)
        for n in range(len(y)):
            dKdz[:,n] = fdiff.second_order_diff(grdata.meank.iloc[:,n], z)
        return dKdy, dKdz
    def calc_mom_transport():
        ddy_upvp = np.zeros(grdims)
        ddz_upwp = np.zeros(grdims)
        d2Udy2 = np.zeros(grdims)
        d2Udz2 = np.zeros(grdims)
        dUdy = np.zeros(grdims)
        dUdz = np.zeros(grdims)
        for n in range(len(z)):
            ddy_upvp[n,:] \
                    = fdiff.second_order_diff(grdata.meanupvp.iloc[n,:], y)
            dUdy[n,:] = fdiff.second_order_diff(grdata.meanu.iloc[n,:], y)
            d2Udy2[n,:] = fdiff.second_order_diff(dUdy[n,:], y)
        for n in range(len(y)):
            ddz_upwp[:,n] \
                    = fdiff.second_order_diff(grdata.meanupwp.iloc[:,n], z)
            dUdz[:,n] = fdiff.second_order_diff(grdata.meanu.iloc[:,n], z)
            d2Udz2[:,n] = fdiff.second_order_diff(dUdz[:,n], z)
        return {"dUdy" : dUdy, "ddy_upvp" : ddy_upvp, "d2Udy2" : d2Udy2,
                "dUdz" : dUdz, "ddz_upwp" : ddz_upwp, "d2Udz2" : d2Udz2}
    if "meanucont" in plotlist or "all" in plotlist:
        # Plot contours of mean streamwise velocity
        plt.figure(figsize=figsize_horiz_contour)
        cs = plt.contourf(y_R, z_H, meanu, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend="both",
                          orientation="horizontal", pad=horiz_cbar_pad)
        cb.set_label(r"$U/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/meanucont"+savetype)
    if "v-wquiver" in plotlist or "all" in plotlist:
        # Make quiver plot of v and w velocities
        plt.figure(figsize=(10,5))
        Q = plt.quiver(y_R, z_H, meanv, meanw)
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
    if "meanu_2tsrs" in plotlist or "all" in plotlist:
        # Plot mean velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanu[ind], marker="o", label=r"$\lambda = 1.9$")
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanu[ind], marker="^", label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$U/U_\infty$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath + "/meanu_2tsrs" + savetype)
    if "meanv_2tsrs" in plotlist or "all" in plotlist:
        # Plot mean cross stream velocity profiles at two TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanv[ind], marker="o", label=r"$\lambda = 1.9$")
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanv[ind], marker="^", label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$V/U_\infty$")
        plt.ylim(-0.1501, 0.1001)
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/meanv_2tsrs"+savetype)
    if "meanw_2tsrs" in plotlist or "all" in plotlist:
        # Plot mean vertical velocity profiles at two TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanw[ind], marker="o", label=r"$\lambda = 1.9$")
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanw[ind], marker="^", label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$W/U_\infty$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath + "/meanw_2tsrs" + savetype)
    if "meanvel_vs_tsr" in plotlist or "all" in plotlist:
        # Plot mean velocity components vs TSR
        tsr = df.tsr.values
        ind1 = [run-1 for run in range(1,32)]
        ind2 = [run-1 for run in range(347,378)]
        plt.figure()
        cycler = plt.rcParams["axes.prop_cycle"]
        markers = ["o", "s", "^"]
        labels = ["$U$", "$V$", "$W$"]
        quantities = [df.meanu, df.meanv, df.meanw]
        for c, m, lab, q in zip(cycler, markers, labels, quantities):
            c = c["color"]
            plt.plot(tsr[ind1], q[ind1], marker=m, markerfacecolor="none",
                     label="", markeredgecolor=c)
            plt.plot(tsr[ind2], q[ind2], marker=m, label=lab, color=c)
        plt.legend(ncol=3)
        plt.xlabel(r"$\lambda$")
        plt.ylabel("Normalized mean velocity")
        styleplot()
        if save:
            plt.savefig(savepath + "/mean_vel_vs_tsr" + savetype)
    if "k_2tsrs" in plotlist or "all" in plotlist:
        # Plot mean velocities at two different TSRs
        plt.figure()
        runs = getruns(z_H=0.25, tsr=1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.k[ind], marker="o", label=r"$\lambda = 1.9$")
        runs = getruns(z_H=0.25, tsr=1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.k[ind], marker="^", label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$k/U_\infty^2$")
        plt.legend(loc="upper right")
        styleplot()
        if save:
            plt.savefig(savepath + "/k_2tsrs" + savetype)
    if "stdu_2tsrs" in plotlist or "all" in plotlist:
        # Plot stdu velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.stdu[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.stdu[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\sigma_u/U_\infty$")
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+"/stdu_2tsrs"+savetype)
    if "uw_2tsrs" in plotlist or "all" in plotlist:
        # Plot uw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanupwp[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanupwp[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{u'w'}$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/uw_2tsrs"+savetype)
    if "stducont" in plotlist or "all" in plotlist:
        # Plot contours of streamwise turbulence intensity
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, grdata["stdu"], 20)
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
    if "uvcont" in plotlist or "all" in plotlist:
        # Plot contours of uv Reynolds stress
        plt.figure(figsize=figsize_horiz_contour)
        cs2 = plt.contourf(y_R, z_H, uv, 15, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=horiz_cbar_pad)
        cb2.set_label(r"$\overline{u^\prime v^\prime}/U_\infty^2$")
        cb2.set_ticks(np.arange(-0.02, 0.025, 0.005), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/uvcont"+savetype)
    if "stdv_2tsrs" in plotlist or "all" in plotlist:
        # Plot stdv velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.stdv[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.stdv[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\sigma_v/U_\infty$")
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+"/stdv_2tsrs"+savetype)
    if "stdw_2tsrs" in plotlist or "all" in plotlist:
        # Plot stdw velocities at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.stdw[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.stdw[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\sigma_w/U_\infty$")
        plt.legend(loc=1)
        styleplot()
        if save:
            plt.savefig(savepath+"/stdw_2tsrs"+savetype)
    if "uv_2tsrs" in plotlist or "all" in plotlist:
        # Plot uv Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanupvp[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanupvp[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{u'v'}$")
        plt.legend()
        styleplot()
        if save:
            plt.savefig(savepath+"/uv_2tsrs"+savetype)
    if "kcont" in plotlist or "all" in plotlist:
        # Plot contours of k
        plt.figure(figsize=figsize_vertical_contour)
        cs = plt.contourf(y_R, z_H, grdata["k"]/(1.0**2), 20,
                          cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1.0, extend="both",
                          orientation="vertical", pad=vertical_cbar_pad)
        cb.set_label(r"$k/U_\infty^2$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/kcont"+savetype)
    if "meankcont" in plotlist or "all" in plotlist:
        # Plot contours of k
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, grdata["meank"]/(0.5*1**2), 20,
                          cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cb.set_label(r"$K/\frac{1}{2}U_\infty^2$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/meankcont"+savetype)
    if "meanvcont" in plotlist or "all" in plotlist:
        # Plot contours of meanv
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, meanv, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cbmv.set_label(r"$\overline{v}/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/meanvcont"+savetype)
    if "stdvcont" in plotlist or "all" in plotlist:
        # Plot contours of stdv
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, grdata["stdv"], 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cb.set_label(r"$\sigma_v/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/stdvcont"+savetype)
    if "meanwcont" in plotlist or "all" in plotlist:
        # Plot contours of meanw
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, meanw, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cbmv.set_label(r"$\overline{w}/U_{\infty}$")
        turb_lines()
        if save:
            plt.savefig(savepath+"/meanwcont"+savetype)
    if "stdwcont" in plotlist or "all" in plotlist:
        # Plot contours of stdw
        plt.figure(figsize=(10,5))
        cmv = plt.contourf(y_R, z_H, grdata["stdw"], 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cbmv = plt.colorbar(cmv, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cbmv.set_label(r"$\sigma_w/U_{\infty}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/stdwcont"+savetype)
    if "vw_2tsrs" in plotlist or "all" in plotlist:
        # Plot vw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanvpwp[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanvpwp[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{v'w'}$")
        plt.legend(loc=4)
        styleplot()
        if save:
            plt.savefig(savepath+"/vw_2tsrs"+savetype)
    if "vwcont" in plotlist or "all" in plotlist:
        # Plot contours of vw Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, vw, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=0.3)
        cb2.set_label(r"$\overline{v^\prime w^\prime}/U_\infty^2$")
        cb2.set_ticks(np.linspace(-.008,.006,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/vwcont"+savetype)
    if "uwcont" in plotlist or "all" in plotlist:
        # Plot contours of vw Reynolds stress
        plt.figure(figsize=figsize_horiz_contour)
        cs2 = plt.contourf(y_R, z_H, uw, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=horiz_cbar_pad)
        cb2.set_label(r"$\overline{u^\prime w^\prime}/U_\infty^2$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/uwcont"+savetype)
    if "vvcont" in plotlist or "all" in plotlist:
        # Plot contours of vv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, vv, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cb2.set_label(r"$\overline{v^\prime v^\prime}$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/vvcont"+savetype)
    if "wwcont" in plotlist or "all" in plotlist:
        # Plot contours of vv Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, ww, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cb2.set_label(r"$\overline{w^\prime w^\prime}$")
#        cb2.set_ticks(np.linspace(-.015,.013,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/wwcont"+savetype)
    if "uucont" in plotlist or "all" in plotlist:
        # Plot contours of uu Reynolds stress
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, grdata.meanupup, 20)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=0.18)
        cb2.set_label(r"$\overline{u^\prime u^\prime}$")
        cb2.set_ticks(np.linspace(0,.108,6), update_ticks=True)
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/uucont"+savetype)
    if "vv_2tsrs" in plotlist or "all" in plotlist:
        # Plot vw Re stress at two different TSRs
        plt.figure()
        runs = getruns(0.25, 1.9)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanvpvp[ind], "-ok", markerfacecolor="none",
                 label=r"$\lambda = 1.9$")
        plt.hold(True)
        runs = getruns(0.25, 1.4)
        ind = [run-1 for run in runs]
        plt.plot(y_R, df.meanvpvp[ind], "--^k", markerfacecolor="none",
                 label=r"$\lambda=1.4$")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$\overline{v^\prime v^\prime}$")
        plt.legend()
        styleplot()
        if save:
            plt.savefig(savepath+"/vv_2tsrs"+savetype)
    if "fpeak_u" in plotlist or "all" in plotlist:
        plt.figure(figsize=(10,5))
        cs2 = plt.contourf(y_R, z_H, grdata["fpeak_u"], cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=0.18)
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
        cs2 = plt.contourf(y_R, z_H, grdata["fstrength_u"], 20,
                           cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=0.18)
        cb2.set_label(r"$S_{\max}/\sigma^2_u$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"/fstrength_u"+savetype)
    if "fpeak_v" in plotlist or "all" in plotlist:
        plt.figure(figsize=figsize_horiz_contour)
        cs2 = plt.contourf(y_R, z_H, grdata["fpeak_v"], cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=horiz_cbar_pad)
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
        plt.figure(figsize=figsize_horiz_contour)
        cs2 = plt.contourf(y_R, z_H, grdata["fstrength_v"], 20,
                           cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=horiz_cbar_pad)
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
        cs2 = plt.contourf(y_R, z_H, grdata["fpeak_w"], cmap=plt.cm.coolwarm,
                           levels=np.linspace(0,10,21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=0.18)
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
        cs2 = plt.contourf(y_R, z_H, grdata["fstrength_w"], 20,
                           cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb2 = plt.colorbar(cs2, shrink=1, extend="both",
                           orientation="horizontal", pad=0.18)
        cb2.set_label(r"$S_{\max}/\sigma^2_w$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath+"fstrength_w"+savetype)
    # Plot estimate for production of turbulence kinetic energy
    if "kprod" in plotlist or "all" in plotlist:
        kprod, meandiss = calc_kprod_meandiss()
        plt.figure(figsize=(10,5))
        cs = plt.contourf(y_R, z_H, kprod, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1, extend="both",
                          orientation="horizontal", pad=0.18)
        cb.set_label(r"$-\overline{u_i' u_j'}\frac{\partial U_i}{\partial x_j}$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        styleplot()
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/kprod"+savetype)
    if "meankadv" in plotlist or "all" in plotlist:
        dKdy, dKdz = calc_meankgrad()
        # Make quiver plot of K advection
        plt.figure(figsize=(10,5))
        plt.hlines(0.5, -1, 1, linestyles="solid", colors="r",
                   linewidth=3)
        plt.vlines(-1, -0.2, 0.5, linestyles="solid", colors="r",
                   linewidth=3)
        plt.vlines(1, -0.2, 0.5, linestyles="solid", colors="r",
                   linewidth=3)
        Q = plt.quiver(y_R, z_H, meanv/meanu*dKdy, meanw/meanu*dKdz,
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
    if "Kturbtrans" in plotlist or "all" in plotlist:
        tt, tty, ttz = calc_meankturbtrans()
        plt.figure(figsize=figsize_horiz_contour)
        cs = plt.contourf(y_R, z_H, tty, 20, cmap=plt.cm.coolwarm,
                          levels=np.linspace(-0.08, 0.08, 21))
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        styleplot()
        cb = plt.colorbar(cs, shrink=1, extend="both",
                          orientation="horizontal", pad=horiz_cbar_pad)
        cb.set_label(r"$-\frac{1}{2}\frac{\partial}{\partial x_j}\overline{u_i^\prime u_j^\prime} U_i$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        if save:
            plt.savefig(savepath+"/Kturbtrans"+savetype)

    if "meancontquiv" in plotlist or "all" in plotlist:
        plt.figure(figsize=figsize_vertical_quiver)
        # Add contours of mean velocity
        cs = plt.contourf(y_R, z_H, meanu, 20, cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=vertical_cbar_shrink, extend="both",
                          orientation="vertical", pad=vertical_cbar_pad)
        cb.set_label(r"$U/U_{\infty}$")
        # Make quiver plot of v and w velocities
        q = plt.quiver(y_R, z_H, meanv, meanw, scale=3, width=0.0022,
                       edgecolor="none")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(q, 0.75, 0.08, 0.1, r"$0.1 U_\infty$",
                      labelpos="E", coordinates="figure",
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
            plt.savefig(savepath+"/meancontquiv"+savetype)

    if "xvorticity" in plotlist or "all" in plotlist:
        z = 1.0*z_H
        y = R*y_R
        dWdy = np.zeros(np.shape(meanu))
        dVdz = np.zeros(np.shape(meanu))
        for n in range(len(z)):
            dWdy[n,:] = fdiff.second_order_diff(meanw.iloc[n,:], y)
        for n in range(len(y)):
            dVdz[:,n] = fdiff.second_order_diff(meanv.iloc[:,n], z)
        # Make quiver plot of K advection
        plt.figure(figsize=figsize_horiz_contour)
        cs = plt.contourf(y_R, z_H, dWdy-dVdz, 20, cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1, extend="both",
                          orientation="horizontal", pad=horiz_cbar_pad)
        cb.set_label(r"$\Omega_x$")
        turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        styleplot()
        if save:
            plt.savefig(savepath + "/xvorticity" + savetype)

    if "Kbargraphs" in plotlist or "all" in plotlist:
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
            quantities = [-meanv[yloc][zloc]/meanu[yloc][zloc]*dKdy[i1,i2],
                          -meanw[yloc][zloc]/meanu[yloc][zloc]*dKdz[i1,i2],
                          tty[i1,i2]/meanu[yloc][zloc],
                          ttz[i1,i2]/meanu[yloc][zloc],
                          kprod[yloc][zloc]/meanu[yloc][zloc],
                          meandiss[i1,i2]/meanu[yloc][zloc]]
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

    if "Kbargraph" in plotlist or "all" in plotlist:
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
        plt.figure(figsize=np.array((7, 3))*scale)
        names = [r"$y$-adv.", r"$z$-adv.",
                 r"$y$-turb.",
                 r"$z$-turb.",
                 r"$k$-prod.", "Mean diss."]
        quantities = [average_over_area(-2*meanv/meanu*dKdy/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(-2*meanw/meanu*dKdz/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*tty/meanu/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*ttz/meanu/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*kprod/meanu/(0.5*U**2)/D, y_R, z_H),
                      average_over_area(2*meandiss/meanu/(0.5*U**2)/D, y_R, z_H)]
        ax = plt.gca()
        ax.bar(range(len(names)), quantities, color=barcolor,
               edgecolor="black", width=0.5)
        ax.set_xticks(np.arange(len(names))+0.25)
        ax.set_xticklabels(names)
        plt.hlines(0, 0, len(names), color="black", linewidth=1)
        plt.ylabel(r"$\frac{K \, \mathrm{ transport}}{UK_\infty D^{-1}}$")
        plt.tight_layout()
        if print_analysis:
            print("K recovery rate (%/D) =",
                  2*np.sum(quantities)/(0.5*U**2)/D*100)
        if save:
            plt.savefig(savepath+"/Kbargraph"+savetype)

    if "mombargraph" in plotlist or "all" in plotlist:
        """Make a bar graph of terms contributing to dU/dx:
          * Cross-stream advection
          * Vertical advection
          * Cross-stream Re stress gradient
          * Vertical Re stress gradient
          * Cross-steam diffusion
          * Vertical diffusion
        """
        data = calc_mom_transport()
        dUdy = data["dUdy"]
        dUdz = data["dUdz"]
        tty = data["ddy_upvp"]
        ttz = data["ddz_upwp"]
        d2Udy2 = data["d2Udy2"]
        d2Udz2 = data["d2Udz2"]
        plt.figure(figsize=np.array((7, 3))*scale)
        names = [r"$-V \frac{\partial U}{\partial y}$",
                 r"$-W \frac{\partial U}{\partial z}$",
                 r"$-\frac{\partial}{\partial y} \overline{u^\prime v^\prime}$",
                 r"$-\frac{\partial}{\partial z} \overline{u^\prime w^\prime}$",
                 r"$\nu \frac{\partial^2 U}{\partial y^2}$",
                 r"$\nu \frac{\partial^2 U}{\partial z^2}$"]
        quantities = [average_over_area(-2*meanv*dUdy/meanu/D, y_R, z_H),
                      average_over_area(-2*meanw*dUdz/meanu/D, y_R, z_H),
                      average_over_area(-2*tty/meanu/D, y_R, z_H),
                      average_over_area(-2*ttz/meanu/D, y_R, z_H),
                      average_over_area(2*nu*d2Udy2/meanu/D, y_R, z_H),
                      average_over_area(2*nu*d2Udz2/meanu/D, y_R, z_H)]
        ax = plt.gca()
        ax.bar(range(len(names)), quantities, color=barcolor, width=0.5,
               edgecolor="black")
        ax.set_xticks(np.arange(len(names))+0.25)
        ax.set_xticklabels(names)
        plt.hlines(0, 0, len(names), color="black", linewidth=1)
        plt.ylabel(r"$\frac{U \, \mathrm{ transport}}{UU_\infty D^{-1}}$")
        plt.tight_layout()
        if print_analysis:
            print("U recovery rate (%/D) =",
                  2*np.sum(quantities)/U/D*100)
        if save:
            plt.savefig(savepath+"/mombargraph"+savetype)


def plot_torque_ripple():
    i = range(31)
    torque_ripple = np.load("Data/Processed/torque_ripple.npy")
    tsr = np.load("Data/Processed/tsr.npy")
    plt.plot(tsr[i], torque_ripple[i], "-ok", markerfacecolor = "none")
    plt.xlabel(r"$\lambda$", labelpad=20)
    plt.ylabel(r"Torque ripple")
    styleplot()
    print("Torque ripple at TSR =", str(tsr[12])+":", torque_ripple[12])


def plot_Re_c():
    a = np.load("Data/Processed/a.npy")
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


def plot_cp(ax=None, fig=None, save=False, savetype=".pdf", savedir="Figures",
            **kwargs):
    """Plot mean power coefficient versus mean tip speed ratio."""
    df = pd.read_csv("Data/Processed/processed.csv").iloc[:31]
    if ax is None:
        fig, ax = plt.subplots()
    if not "marker" in kwargs.keys():
        kwargs["marker"] = "o"
    ax.plot(df.tsr, df.cp, **kwargs)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$C_P$")
    if fig is not None:
        fig.tight_layout()
    if save:
        if not os.path.isdir(savedir):
            os.makedirs(savedir)
        fig.savefig(os.path.join(savedir, "cp" + savetype))


def plot_cd(ax=None, fig=None, save=False, savetype=".pdf", savedir="Figures",
            **kwargs):
    """Plot mean drag coefficient versus mean tip speed ratio."""
    df = pd.read_csv("Data/Processed/processed.csv").iloc[:31]
    if ax is None:
        fig, ax = plt.subplots()
    if not "marker" in kwargs.keys():
        kwargs["marker"] = "o"
    ax.plot(df.tsr, df.cd, **kwargs)
    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$C_D$")
    ax.set_ylim((0, 1.2))
    if fig is not None:
        fig.tight_layout()
    if save:
        if not os.path.isdir(savedir):
            os.makedirs(savedir)
        fig.savefig(os.path.join(savedir, "cd" + savetype))


def plotperf(plotlist=["cp", "cd"], subplots=True, save=False,
             savepath="Figures", savetype=".pdf", print_perf=True, **kwargs):
    i = np.arange(31)
    data = pd.read_csv("Data/Processed/processed.csv")
    cp = data.cp
    cd = data.cd
    tsr = data.tsr
    eta2 = data.eta2
    if not subplots:
        # Exergy efficiency plot
        plt.figure()
        plt.plot(tsr[i], eta2[i], marker="o", **kwargs)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$\eta_{II}$")
        styleplot()
        if save:
            plt.savefig(savepath+"/eta2"+savetype)
        # Power coefficient plot
        plt.figure()
        plt.plot(tsr[i], cp[i], marker="o", **kwargs)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
        styleplot()
        if save:
            plt.savefig(savepath + "/cpvstsr" + savetype)
        # Drag coefficient plot
        plt.figure()
        plt.plot(tsr[i], cd[i], marker="o", **kwargs)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_D$")
        plt.ylim(0, 1.2001)
        styleplot()
        if save:
            plt.savefig(savepath + "/cdvstsr" + savetype)
        # Torque coefficient plot
        ct = cp/tsr
        plt.figure()
        plt.plot(tsr[i], ct[i], "-ok", **kwargs)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_T$")
        styleplot()
        if save:
            plt.savefig(savepath+"/ctvstsr"+savetype)
    else:
        plt.figure(figsize=(7.5, 3.25))
        plt.subplot(121)
        plt.plot(tsr[i], cp[i], marker="o", **kwargs)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
        plt.grid(True)
        plt.subplot(122)
        plt.plot(tsr[i], cd[i], marker="o", **kwargs)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_D$")
        plt.ylim(0, 1.2001)
        styleplot()
        if save:
            plt.savefig(savepath + "/perf" + savetype)
    if print_perf:
        print("At tsr = 1.9, C_P =",
              cp[np.where(np.round(tsr, decimals=2)==1.9)[0][0]],
              "; C_D =", cd[np.where(np.round(tsr, decimals=2)==1.9)[0][0]])


def plotperf_periodic():
    i = range(31)
    d = pd.read_csv("Data/Processed/processed.csv")
    plt.figure()
    plt.plot(d.tsr[i], d.amp_cd[i])
    styleplot()
    plt.figure()
    plt.plot(d.tsr[i], d.phase_cd[i])


def plot_phase_average(run=13, plot_cp=True, plot_cd=False):
    t1 = 13
    t2 = 30
    t, angle, Ttrans, Tarm, drag, rpm, tsr = loadtdms(run)
    omega = rpm*2*np.pi/60
    cp = Ttrans*omega/(1/2*rho*A_t*U**3)
    cd = drag/(1/2*rho*A_t*U**2)
    angle1 = angle[t1*2000]
    undershoot = 360-np.mod(angle1,360)
    angle1 = angle1 + undershoot
    angle2 = angle[t2*2000]
    overshoot = np.mod(angle2,360)
    angle2 = angle2 - overshoot
    nrevs = (angle2 - angle1)/360
    def find_index(angle, ta):
        i = np.where(np.round(angle, decimals=0)==ta)
        i = i[0]
        if len(i) > 1: i = i[0]
        return i
    i1 = find_index(angle, angle1)
    i2 = find_index(angle, angle1+360)
    npoints = i2-i1
    cp_phave = cp[i1:i2]
    cd_phave = cd[i1:i2]
    print("Averaging over {} revolutions".format(nrevs))
    for n in range(1,int(nrevs)):
            ta1 = angle1+n*360
            i1 = find_index(angle, ta1)
            i2 = i1+npoints
            cp_phave += cp[i1:i2]
            cd_phave += cd[i1:i2]
    cp_phave /= nrevs
    cd_phave /= nrevs
    angleb = np.linspace(0, 360, num=npoints)
    if plot_cp:
        plt.plot(angleb, cp_phave, "k")
    if plot_cd:
        plt.plot(angleb, cd_phave)
    plt.xlabel(r"$\theta$ (deg)")
    plt.ylabel(r"$C_P$")
    styleplot()


if __name__ == "__main__":
    pass
