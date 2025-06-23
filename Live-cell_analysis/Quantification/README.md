# Quantification
LiveCell_Analysis.py quantifies distance time series of TAD anchors. The datasets can be downloaded on Zenodo.
The script allows to compute the following (in about 2 min):
-	Proximal state segmentation
-	Proximal state fraction
-	Proximal state frequency
-	Proximal state lifetime
-	Closing rate
-	Fractions of proximal, extruding and open states derived from an analytical model of coordinate difference distribution.

Other files contain functions used in LiveCell_Analysis.py.

## Requirements

This package is developed in python and requires the following libraries. You can create a conda environment (less than 1 min) as follows:

* conda create -n sabate_env python=3.11.7
* conda activate sabate_env
* pip install pandas==2.1.4 matplotlib==3.8.0 numpy==1.26.4 

The distance time series datasets can be downloaded in Zenodo.
