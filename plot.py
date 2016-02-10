#!/usr/bin/env python
"""This script generates figures from the RVAT baseline experiment."""

from pyrvatbl.plotting import *
import argparse
import sys


jotfigs = ["meancontquiv", "xvorticity", "fpeak_v", "fstrength_v", "uvcont",
            "uwcont", "Kturbtrans", "kcont", "Kbargraph", "mombargraph"]

wakefigs_all = jotfigs + ["meanu_2tsrs", "meanv_2tsrs", "meanw_2tsrs",
                          "meanvel_vs_tsr", "k_2tsrs"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create figures from the "
                                     "UNH-RVAT baseline experiment.")
    parser.add_argument("figures", nargs="*", help="Which figure to create",
                        choices=wakefigs_all + ["multispec", "perf", []],
                        default=None)
    parser.add_argument("--no-subplots", action="store_true",
                        help="Don't use subplots for performance plots")
    parser.add_argument("--all", "-a", help="Plot all figures",
                        action="store_true", default=False)
    parser.add_argument("--no-errorbars", "-e", help="Do not plot error bars",
                        action="store_true", default=False)
    parser.add_argument("--save", "-s", help="Save figures",
                        action="store_true", default=False)
    parser.add_argument("--savetype", help="Format to save figures",
                        default=".pdf")
    parser.add_argument("--style", help="Matplotlib stylesheet")
    parser.add_argument("--no-show", help="Do not show figures",
                        action="store_true", default=False)
    parser.add_argument("--verbose", "-v", help="Print analysis",
                        action="store_true", default=False)
    args = parser.parse_args()

    if not args.figures and not args.all:
        parser.print_help()
        sys.exit("No figures chosen to plot")

    if args.style:
        plt.style.use(args.style)
    else:
        from pxl.styleplot import set_sns
        set_sns()
    savetype = args.savetype
    save = args.save
    errorbars = not args.no_errorbars
    subplots = not args.no_subplots
    if save:
        if not os.path.isdir("Figures"):
            os.makedirs("Figures")

    if not isinstance(args.figures, list):
        args.figures = [args.figures]
    wakefigs = [f for f in args.figures if f in wakefigs_all]

    if "perf" in args.figures or args.all:
        plotperf(subplots=subplots, save=save, savetype=savetype,
                 print_perf=args.verbose)
    if wakefigs or args.all:
        if args.all:
            wakefigs = wakefigs_all
        plotwake(wakefigs, save=save, print_analysis=args.verbose)
    if "multispec" in args.figures or args.all:
        plotmultispec(n_band_average=4, save=save, plot_conf_int=errorbars)

    if not args.no_show:
        plt.show()
