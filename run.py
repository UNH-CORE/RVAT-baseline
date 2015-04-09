#!/usr/bin/env python
"""
This script generates figures from the RVAT baseline experiment.
"""

from Modules.plotting import *
from Modules.processing import *

setpltparams(seaborn=True, latex=False)
save = True
show = True

def main():    
    # List of plots for JoT paper
    jotplots = ["meancomboquiv", "xvorticity", "fpeak_v", "fstrength_v",
                "uvcont", "uwcont", "Kturbtrans", "kcont", "Kbargraph",
                "mombargraph"]
    
    plotperf(subplots=True, save=save)
    plotwake(jotplots, save=save, print_analysis=True)
    plotmultispec(n_band_average=4, save=save, plot_conf_int=True)
    
    if show:
        plt.show()

if __name__ == "__main__":
    if not os.path.isdir("./Figures"):
        os.mkdir("Figures")
    main()
