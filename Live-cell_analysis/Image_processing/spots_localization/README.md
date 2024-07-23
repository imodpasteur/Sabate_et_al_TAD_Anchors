# Spots localization

This repository contains four files for chromatic aberration correction and spot localization on 2-colors images

## License

MIT

The plugin is distributed under the license MIT on an "as is" basis, without warranties or conditions of any kind.

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


1. replicatedSpotDetection.py takes a tif image as input, and outputs a tif file with localizations. The intensities of the pixels in the output image indicate whether the localized spots are elongated or not.

2. shift_correction.py takes a tif image as input, and outputs a tif image at the same location corrected for temporal shifts. The correction is pixelic. No interpolation is performed. This correction allows to better track spots.

3. Use Trackmate Fiji plugin to detect and track red and green spots.

4. trackmate_color_registration.py corrects the chromatic aberrations from the localizations. It uses as input the track file given by Trackmate and the chromatic aberration parameters computed using Chromagnon.

5. relocalization.py allows to improve Trackmate localization by fitting a Gaussian model to spots, and estimates the Cramér Rao bounds, that provide us the localization precision for each localization. It takes as input a tif file, the tracks given by trackmate, and chromatic aberration parameters computed using Chromagnon.
