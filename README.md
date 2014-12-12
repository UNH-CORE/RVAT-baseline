# UNH-RVAT baseline performance and near-wake measurements

This repository contains the processing and plotting code, as well as the 
derived data set from the UNH-RVAT tow tank measurements performed in Spring 2013.

Download/usage
--------------

To download, use the git clone URL to the right. For example, in a terminal

    git clone https://github.com/UNH-CORE/UNH-RVAT-baseline.git

To generate plots, use functions in the `plotting` module. See this
[IPython notebook](http://nbviewer.ipython.org/github/UNH-CORE/UNH-RVAT-baseline/blob/master/Documents/examples.ipynb "View on nbviewer.ipython.org") for examples.

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
@Misc{Bachant2014_data,
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

### Commit `f8cc166`
```bibtex
@INPROCEEDINGS{Bachant2013
  author = {Bachant, P. and Wosnik, M.},
  title = {Performance and near-wake measurements for a vertical axis turbine
	   at moderate Reynolds number},
  booktitle = {Proceedings of the ASME 2013 Fluids Engineering Division Summer Meeting},
  year = {2013},
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
