# Spots localization

This repository contains five files for chromatic aberration correction and spot localization on 2-colors images


## Requirements

This package is developed in python and requires the following libraries

numpy
pandas
scipy
scikit-learn
scikit-image
tifffile
statsmodels
pathlib
imageio
pillow


## Usage

The files are executed in the following order in order to track spots in 2-color imaging. The images have to be in tif format in the following order: TZCYX.


1. replicatedSpotDetection.py takes a 5D tif image as input, and outputs a 2D histogram of localizations saved in tif format. The intensities of the pixels in the output image indicate whether the localized spots are elongated or not. This script was used to filter out replicated spots that look like elongated spots on images.

2. shift_correction.py takes a tif image as input, and outputs a tif image at the same location corrected for temporal shifts. The correction is pixelic. No interpolation is performed. This correction allows to better track spots.

3. Use [Trackmate](https://github.com/trackmate-sc) [Fiji](https://imagej.net/software/fiji/) plugin to detect and track red and green spots.

4. trackmate_color_registration.py is then used to correct chromatic aberrations. It takes as input registration parameters provided by [Chromagnon](https://github.com/macronucleus/Chromagnon) and the xml file provided by Trackmate.

5. Run the color pairing using ??? Fiji plugin.

5. Finally, relocalization.py allows to improve trackmate localization by fitting a Gaussian model to spots, and estimates the Cramér Rao bounds, that provide us the localization precision for each localization. It takes as input a tif file (with shifts corrected), the xml file provided by trackmate and the registration parameters provided by Chromagnon.
