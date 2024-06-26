# xas_normalisation
Scripts for normalising and applying monochromator theta offsets for BM31 XAS data. (SNBL, ESRF).

Use: clone or download the repository. Run 'pip install -e .' inside. 

This creates a program called 'xasCorn' (xas COlumn selection Regridding and Normalisation) which can be run from the directory of choice. First it separates the scans from a file into separate files in a subdirectory called 'columns' by finding the ion chambers with the highest counts, then checking if they and the fluorescence counts are above a certain value. Then it regrids them and puts the regridded files into a 'regrid' subdirectory, then it normalises the regridded files into 'norm' and 'merge' folders. After 1 cycle it enters a loop looking for new files or if files have been modified so it can run while measuring. Use -to/--thetaoffset to apply an angular correction to the monochromator, e.g. when using the side crystal, and --unit to choose whether to put the regrid and normalisation files in eV or keV.

xasCorn --help<br>
usage: xasCorn [-h] [-to THETAOFFSET] [-u UNIT] [directory]

positional arguments:
  directory             the directory to run in. Default is current directory

options:
  -h, --help            show this help message and exit<br>
  -to THETAOFFSET, --thetaOffset THETAOFFSET
                        theta offset to apply to monochromator for energy correction<br>
  -u UNIT, --unit UNIT  keV or eV (default keV) for the regrid and normalised files

It uses Larch (https://xraypy.github.io/xraylarch/index.html, https://iopscience.iop.org/article/10.1088/1742-6596/430/1/012007) to normalise the data.
