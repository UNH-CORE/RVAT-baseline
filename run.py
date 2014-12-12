#!/usr/bin/env python
"""
This script is used to call processing and plotting functions.

"""

from Modules.plotting import *


def main():
    setpltparams(latex=True)
    plt.close("all")
    p = "Google Drive/Research/Papers/JoT CFT near-wake/Figures"
    if "linux" in sys.platform:
        p = "/home/pete/" + p
    elif "win" in sys.platform:
        p = "C:/Users/Pete/" + p
        
    jotplots = ["meancomboquiv", "xvorticity", "fpeak_v", "fstrength_v",
                "uvcont", "uwcont", "Kturbtrans", "kcont", "Kbargraph",
                "mombargraph"]
        
#    plotsinglerun(41, perf=True, wake=False, autocorr=False, xaxis="angle")
#    plot_phase_average(124)
#    plotvelspec(y_R=1.5, z_H=0.25, tsr=1.9, show=True)
#    plotperfspec(y_R=1.5, z_H=0.25, tsr=1.9, show=True)
#    plotperf(subplots=True, save=True, savepath=p)
    plotwake(["Kbargraph", "mombargraph"], save=False, savepath=p, print_analysis=True)
#    plotmultispec(n_band_average=5, save=True, savepath=p)
#    plotperf_periodic()
#    plotvelhist(5)

if __name__ == "__main__":
    main()