# xas_normalisation
Scripts for normalising and applying monochromator theta offsets for BM31 XAS data. (SNBL, ESRF).

Use: clone or download the repository. Run 'pip install -e .' inside. 

This creates a program called 'xasCorn' (xas COlumn selection Regridding and Normalisation) which can be run from the directory of choice. First it separates the scans from a file into separate files in a subdirectory called 'columns' by finding the ion chambers with the highest counts, then checking if they and the fluorescence counts are above a certain value. Then it regrids them and puts the regridded files into a 'regrid' subdirectory, then it normalises the regridded files into 'norm' and 'merge' folders. After 1 cycle it enters a loop looking for new files or if files have been modified so it can run while measuring. Use -to/--thetaoffset to apply an angular correction to the monochromator, e.g. when using the side crystal, and --unit to choose whether to put the regrid and normalisation files in eV or keV.
```
usage: xasCorn [-h] [-to THETAOFFSET] [-u UNIT] [-av AVERAGING] [-e ELEMENTS] [-ee EXCLUDEELEMENTS] [-s SUBDIR]
               [-d DSPACING]
               [directory]

positional arguments:
  directory             the directory to run in. Default is current directory

options:
  -h, --help            show this help message and exit
  -to, --thetaOffset THETAOFFSET
                        theta offset to apply to monochromator for energy correction. For side crystal approximately
                        -0.103 with default d-spacing 3.13429
  -u, --unit UNIT       keV or eV (default keV) for the regrid and normalised files
  -av, --averaging AVERAGING
                        number of spectra to average over, default 1 (no averaging)
  -e, --elements ELEMENTS
                        comma separated list of elements, e.g. "-e Fe,Cu", if nothing given will search all elements
  -ee, --excludeElements EXCLUDEELEMENTS
                        comma separated list of elements to exclude, e.g. "-ee Fe,Cu", if nothing given will not
                        exculde any (doesn't do anything if used with -e)
  -s, --subdir SUBDIR   how to arrange output folders. "edge" - everything with the same element and type (EXAFS or
                        XANES) will go in the same folder, or "file" - every file gets its own folder. Default "edge"
  -d, --dspacing DSPACING
                        monochromator d-spacing (default 3.13429). If data from before 9/2025, should be 3.13379 (or
                        use the new value and apply a theta offset, try -0.005198). Recalibrated to 3.13429 9/2025, so
                        use this after 8/2025
```
It uses Larch (https://xraypy.github.io/xraylarch/index.html, https://iopscience.iop.org/article/10.1088/1742-6596/430/1/012007) to normalise the data.
