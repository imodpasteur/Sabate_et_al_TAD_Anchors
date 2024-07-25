# Spot localization

This repository contains files for chromatic aberration correction and spot localization on 2-channel images


## Requirements

This package is developed in python and requires the following libraries:

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

The files are executed in the following order to track spots in 2-channel images. The images have to be in tif format in the following order: TZCYX.


1. replicatedSpotDetection.py takes a 5D .tif image as input and outputs a 2D histogram of localizations saved in .tif. The intensities of the pixels in the output image indicate whether the localized spots are elongated or not. This script was used to guide the filtering of replicated spots.

2. shift_correction.py takes a tif image as input, and outputs a .tif image at the same location corrected for temporal XY shifts. The correction is at the pixel resolution. No interpolation is performed. This correction allows to optimize the tracking of spots.

3. Use [Trackmate](https://github.com/trackmate-sc) [Fiji](https://imagej.net/software/fiji/) plugin to detect and track red and green spots independently.

4. trackmate_color_registration.py is then used to correct the coordinates of fluorescent spots for chromatic aberrations. It takes as input registration parameters provided by [Chromagnon](https://github.com/macronucleus/Chromagnon) and the .xml tracking file from Trackmate.

5. Run the 'Pair TrackMate files' Fiji plugin to pair tracks from each channel together, using the tracks corrected for chromatic aberrations.

5. Finally, relocalization.py allows to improve trackmate localization by fitting a Gaussian model to fluorescent spots, and estimates the Cramér Rao bounds. These bounds define the localization precision of each localization. It takes as input a .tif file (with shifts correction), the .csv file from the 'Pair TrackMate files' plugin and the registration parameters provided by Chromagnon. Before relocalization of spots, we linearly interpolated the positions of missing spots, using the known positions of spots before and after missing spots.
