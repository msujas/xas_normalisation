# xas_normalisation
Scripts for normalising and applying monochromator theta offsets for BM31 XAS data. (SNBL, ESRF).

First run the columnExtraction_2thetaCorrection.py script, inputting the root directory at the top. If there's an offset on the monochromator theta angle
(this will likely be the case when using the side crystal), an offset can be set, and it will correct the energy axis. After that it does the job of the 
Column Selector program that is also used on the beamline - separating out scans from individual files and putting them in one sub folder for each file. It goes through all 
sub-folders to find xas data (.dat files). It then regrids the data onto a common axis (sometimes a few points are off) creating 'regrid' subfolders.

The xasNormalisation.py script should be run afterwards. It uses Larch (https://xraypy.github.io/xraylarch/index.html) to batch normalise all xas data in the 
'regrid' folders created by the columnExtraction_2thetaCorrection.py script.
