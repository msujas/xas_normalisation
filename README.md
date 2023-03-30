# xas_normalisation
Scripts for normalising and applying monochromator theta offsets for BM31 XAS data. (SNBL, ESRF).

First run the columnExtraction_thetaCorrection.py script, inputting the root directory at the top. If there's an offset on the monochromator theta angle
(this will likely be the case when using the side crystal), an offset can be set, and it will correct the energy axis. After that it does the job of the 
Column Selector program that is also used on the beamline - separating out scans from individual files and putting them in one sub folder for each file. 
It goes through all sub-folders to find xas data (.dat files). It then regrids the data onto a common axis (sometimes a few points are off) creating 
'regrid' subfolders. It uses the BM31 counter names for energy, ion chambers, monochromator theta angle and fluorescence. It only chooses monitor and I1 columns which have a maximum value > 10000 * step time (searching for columns with mon_ and ion_1 in the names).

The xasNormalisation.py script should be run afterwards. It uses Larch (https://xraypy.github.io/xraylarch/index.html, https://iopscience.iop.org/article/10.1088/1742-6596/430/1/012007) to batch normalise all xas data in the 'regrid' folders created by the columnExtraction_thetaCorrection.py script. Monitor and transmission ion chamber, and fluorescence counter names are those used for BM31. It will search for columns with mon_ and ion_1 in the names to determine which is correct.
