

path_image_5D=./example/example.tif
path_chromagnon=./example/example_chromagnon.csv



python replicatedSpotDetection.py $path_image_5D $path_chromagnon





python shift_correction.py $path_image_5D $path_chromagnon
path_drift="./example/example_drift.csv"
path_image_DC="./example/example_DC.tif" #path shift corrected image (NO interpolation)
path_image_DC_CHAB="./example/example_DC_CHAB.tif" #path shift and CHromaticABerrations corrected image an (INCLUDING interpolations : used for visualization, not for localization)




echo "Open example_DC.tif image in Fiji and run Trackmate to detect spots on both Red and green channels and save as example_red.xml and example_green.xml files"
path_red_trackmate="./example/example_red.xml"
path_green_trackmate="./example/example_green.xml"




python trackmate_color_registration.py $path_red_trackmate $path_green_trackmate $path_chromagnon $path_drift $path_image_DC_CHAB
path_red_track_registered="./example/example_red_CHAB.xml"
path_green_track_registered="./example/example_green_CHAB.xml"


echo "Open Fiji and run Pair_Trackmate_files plugin to pair example_red_CHAB.xml and example_green_CHAB.xml files and save as example_PAIRED.csv"
path_paired_tracks="./example/example_PAIRED.csv"



python relocalization.py $path_image_DC $path_chromagnon $path_paired_tracks $path_drift
path_result_tracks="./example/example_PAIRED_GaussianFit.csv"




python ./example/plot_example_result.py $path_result_tracks
