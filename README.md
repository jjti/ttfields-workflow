# TTFields-Workflow
An SPM8 toolbox for automating MRIs to FEA models for TTFields

### Installation and use
1. Download the repo as a zip file and unzip into SPM8/toolbox (or elsewhere)
2. Rename to TTFields-Workflow (removing the branch name from the folder)
3. Within MATLAB, add the SPM8 path if it has not been already; eg: "addpath('E:\spm8')"
4. Add the path and subpaths of TTFields-Workflow; eg: "addpath(genpath('E:\spm8\toolbox\TTfields-Workflow'))"
5. Type and enter "TTFs" in the MATLAB console
6. Within the prompt, navigator to a DICOM directory for segmentation and select all images; ctrl-A then Enter

### Optional function calls
pyScanIPScript and batScript, uncommented, can be used to generate ScanIP post-processing and initialization scripts, respectively

### Dependencies
The dependencies have been included, but they are listed below along with their sources

##### [Huang et al., 2013](http://bme.ccny.cuny.edu/faculty/lparra/autosegment/)
Y. Huang, J.P. Dmochowski, Y. Su, A. Datta, C. Rorden, L.C. Parra, Automated MRI Segmentation for Individualized Modeling of Current Flow in the Human Head, Journal of Neural Engineering, 10(6):066004, 2013.

##### [NIfTI_tools](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
Tools for NIfTI and ANALYZE image by [Jimmy Shen](https://www.mathworks.com/matlabcentral/profile/authors/757722-jimmy-shen)

##### [cylinder2P](https://www.mathworks.com/matlabcentral/fileexchange/21758-cylinder-surface-connecting-2-points?focused=5104454&tab=function)
Cylinder generation function by [Wei Pan](https://www.mathworks.com/matlabcentral/profile/authors/1453399-wei-pan)

##### [knnsearch](https://www.mathworks.com/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit?focused=5151612&tab=function)
k-nearest neighbors search by [Yi Cao](https://www.mathworks.com/matlabcentral/profile/authors/69713-yi-cao)

##### [normnd](http://www.mathworks.com/matlabcentral/fileexchange/41609-point-cloud-normal-vector?focused=3785579&tab=function)
normal vector algo by [Jered Wells](https://www.mathworks.com/matlabcentral/profile/authors/2714347-jered-wells)
