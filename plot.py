#!/usr/bin/env python
"""
This script generates figures from the RVAT baseline experiment.
"""

from py_rvat_baseline.plotting import *

setpltparams(seaborn=True, latex=False)
save = True
show = True

def main():
    # List of plots for JoT paper
    jotplots = ["meancontquiv", "xvorticity", "fpeak_v", "fstrength_v",
                "uvcont", "uwcont", "Kturbtrans", "kcont", "Kbargraph",
                "mombargraph"]
    plotperf(subplots=False, save=save)
    plotperf(subplots=True, save=save)
    plotwake(jotplots, save=save, print_analysis=True)
    plotmultispec(n_band_average=4, save=save, plot_conf_int=True)

    if show:
        plt.show()

if __name__ == "__main__":
    if not os.path.isdir("./Figures"):
        os.mkdir("Figures")
    main()
