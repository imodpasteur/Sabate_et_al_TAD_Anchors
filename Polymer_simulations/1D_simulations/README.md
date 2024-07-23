# Spots localization

This repository allows to simulate in 1D a topologically associating domain (TAD) of a given size. It simulates cohesin mediated loops randomly loaded along the chromatin with a given density. Loops extrude chromatin at a given velocity for a given processivity in average. The loop can stop at CTCF locations for a time defined in CTCF_positions folder and determined using CHIP-seq data.

## License

MIT

These scripts are distributed under the license MIT on an "as is" basis, without warranties or conditions of any kind.

## Requirements

This package is developed in python and requires the following libraries

numpy Pillow scipy pandas ipython scikit-image scikit-learn imageio tifffile pycairo seaborn igraph


## Usage

The script runSimulationBeforeLammps.py allows to reproduce 1D simulations as it was done in the paper. This script uses as input:

* A processivity (we varied it from 50 to 2000 Kb in our paper)
* A density of cohesins (we varied it from 2 to 40 (coh/Mb) in our paper)
* An averaged CTCF attach time (we set 150 s in our paper)
* The size of the cell line among the four used in our paper {345, 566, 918, 576 (Kb)}
* An extrusion velocity (we varied it from 0.125 to 0.5 Kb/s in our paper)
* The propability of CTCF occupancy (we used 0.5 in our paper)
* The output simulation path where the simulations are stored

The simulations were then used in Lammps to simulate Polymer.



To simulate new 1D polymers with your own set of parameters, you may use the class Chromatin.py. It allows to simulate a piece of chromatin in 1D. Run example.py file in order see how to use it.
The Chromatin.py constructor uses as input :
* processivityKb (The average amount of extrusion per loop)
* cohesinNumberPerMb (the density of loops)
* extrusionSpeedKb=1 (the extrusion velocity)
* avg_CTCF_attach_time (The time CTCF sites are attached to chromatin)
* proba_CTCF_occupancy (The probability that a CTCF is present when a loop meet its position)
* deltaT (The time step)
* totalAreaMb (the size of chromatin)
The following methods allows then to simulate loop extrusion
* addCTCFsite() allows to add a CTCF site at a position with a given probability to block cohesin
* setAnchors() allows to set the position of two anchors. The simulation will compute the distance between these anchors at each step (deltaT)
* run() allows to run the simulation for a given number of steps
* computeAnchor2AnchorDistance() allows to compute the distance between the two anchors
* computeAndGet_HiC_map() allows to compute the HiC map
