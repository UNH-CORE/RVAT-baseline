# UNH-RVAT baseline performance and near-wake measurements

This repository contains the processing and plotting code, as well as the 
derived data set from the UNH-RVAT tow tank measurements performed in Spring 2013.

Download/usage
--------------

To download, use the git clone URL to the right. For example, in a terminal

    git clone https://github.com/UNH-CORE/UNH-RVAT-baseline.git

To generate plots, execute `run.py` with Python or IPython. Also see this
[IPython notebook](http://nbviewer.ipython.org/github/UNH-CORE/RVAT-baseline/blob/master/Documents/examples.ipynb "View on nbviewer.ipython.org") for examples.

To contribute back to the main repository, use GitHub's fork/pull mechanism.

### Dependencies

  * NumPy
  * SciPy
  * matplotlib
  * pandas
  * [pyTDMS](https://github.com/petebachant/pytdms) (for processing raw data)
  * [PXL](https://github.com/petebachant/PXL)

To install all:

```
pip install -r requirements.txt
```

## How to cite
Please cite 

```bibtex
@Misc{Bachant2014-RVAT-baseline,
  Title                    = {UNH-RVAT baseline performance and near-wake measurements: Reduced dataset and processing code},
  Author                   = {Peter Bachant and Martin Wosnik},
  HowPublished             = {fig\textbf{share}. http://dx.doi.org/10.6084/m9.figshare.1080781},
  Month                    = {June},
  Year                     = {2014},
  Doi                      = {10.6084/m9.figshare.1080781},
  Url                      = {http://dx.doi.org/10.6084/m9.figshare.1080781}
}
```

Publications
------------
These data were used in the following publications:

### Commit `478fdd5`

Bachant, P. and Wosnik, M. (2015) 
[Characterising the near-wake of a cross-flow turbine]
(http://doi.org/10.1080/14685248.2014.1001852). _Journal of Turbulence_

```bibtex
@Article{Bachant2015-JoT,
  Title                    = {Characterising the near-wake of a cross-flow turbine},
  Author                   = {Peter Bachant and Martin Wosnik},
  Journal                  = {Journal of Turbulence},
  Year                     = {2015},
  Month                    = {January},
  Number                   = {4},
  Pages                    = {392-410},
  Volume                   = {16},
  Doi                      = {10.1080/14685248.2014.1001852}
}
```

### Commit `f8cc166`

Bachant, P. and Wosnik, M. (2013) [Performance and near-wake measurements
for a vertical axis turbine at moderate Reynolds number]
(http://doi.org/10.1115/FEDSM2013-16575).

```bibtex
@InProceedings{Bachant2013,
  Title                    = {Performance and near-wake measurements for a vertical axis turbine at moderate {R}eynolds number},
  Author                   = {Bachant, Peter and Wosnik, Martin},
  Booktitle                = {Proceedings of the ASME Fluids Engineering Division Summer Meeting},
  Year                     = {2013},
  Address                  = {Incline Village, NV},
  Month                    = {July},
  Number                   = {FEDSM2013-16575},
  Doi                      = {10.1115/FEDSM2013-16575}
}
```

Other resources
---------------

Turbine CAD (STEP) files are available at http://figshare.com/articles/UNH_RVAT_CAD_models/1062009

License
-------
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/4.0/88x31.png" />
</a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
Creative Commons Attribution 4.0 International License</a>.
