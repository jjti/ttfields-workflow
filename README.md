# TTFields-Workflow
An SPM8 toolbox for automating MRIs to FEA models for TTFields

### Installation and use
Download the repo as a zip file and unzip into SPM8/toolbox
Add the SPM8 path (addpath)
TTFs is the primary function (triggers a prompt for selecting the DICOM directory to segment)

### Optional function calls
pyScanIPScript and batScript, uncommented, can be used to generate ScanIP post-processing and initialization scripts, respectively

### Dependencies
Several dependencies are needed, and they have been included, but they are listed below along with their sources:

##### [Huang et al., 2013](http://bme.ccny.cuny.edu/faculty/lparra/autosegment/)
Y. Huang, J.P. Dmochowski, Y. Su, A. Datta, C. Rorden, L.C. Parra, Automated MRI Segmentation for Individualized Modeling of Current Flow in the Human Head, Journal of Neural Engineering, 10(6):066004, 2013.

##### [NIfTI_tools](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
MATLAB Central
Tools for NIfTI and ANALYZE image by Jimmy Shen

##### [cylinder2P](https://www.mathworks.com/matlabcentral/fileexchange/21758-cylinder-surface-connecting-2-points?focused=5104454&tab=function)
MATLAB Central
by Wei Pan

##### [knnsearch](https://www.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit?focused=5151612&tab=function)
MATLAB Central
by Yi Cao

##### [normnd](http://www.mathworks.com/matlabcentral/fileexchange/41609-point-cloud-normal-vector?focused=3785579&tab=function)
MATLAB Central
by Jered Wells

