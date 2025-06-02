# 1D simulations

This repository allows to simulate loop extrusion on a genomic region of a given size in 1D. It simulates cohesin-mediated loops randomly loaded along the polymer with a given density, velocity and average processivity. The 'CTCF_positions' folder defines the location, orientation and normalized median ChIP-Seq fold-enrichment of CTCF sites along the simulated polymer, together with the location of TAD anchors.


## Requirements

This package is developed in python and requires the following libraries. You can create a conda environment as follow to run the program:

* conda create -n sabate_repo python=3.10
* conda activate sabate_repo
* pip install imageio numpy scipy scikit-learn tifffile elementpath scikit-image pandas matplotlib seaborn pycairo python-igraph statsmodels


## Usage

The script runSimulationBeforeLammps.py allows to reproduce 1D simulations as it was done in the paper. This script uses as input:

* A cohesin processivity (varied from 50 to 2000 kb in our paper). The processivity is set on naked DNA. The effective genomic processivity is lower and depends on the genomic region due to stalling at CTCF sites
* A cohesin density (varied from 2 to 40 cohesin/Mb in our paper)
* An average CTCF residence time (set to 150 s in our paper). This represents the median CTCF residence time genome-wide. The normalized ChIP-Seq fold-enrichment is multiplied by this aveerage to obtain the CTCF residence time associated with each CTCF site.
* The size of the TAD (345, 566, 918, 576 kb in our paper)
* An extrusion motor speed (we varied it from 0.0625 to 0.5 kb/s in our paper). This is the speed of one-directional extrusion. The effective speed is therefore twice the value indicated due to bidirectional extrusion.
* The propability of CTCF occupancy (we used 0.5 in our paper)
* The output simulation path where the simulations are stored
These 1D simulations were then used as inputs for LAMMPS to simulate 3D polymer dynamics, together with loop extrusion.
As an example, we provide "runSimulationAsInPaper.sh" file that show how to run simulations. It creates a set of folders (one folder for each condition) including "run_Multiple_loop.sh" that can then be run on a cluster for polymer simulation with LAMMPS.




To simulate new 1D polymers with your own set of parameters, you may use the class Chromatin.py. It allows to simulate a piece of chromatin in 1D. Run example.py file in order to see how to use it.
The Chromatin.py constructor uses as input:
* processivityKb (average amount of extrusion per loop, in kb)
* cohesinNumberPerMb (density of cohesin per Mb)
* extrusionSpeedKb (unidirectional cohesin motor speed, in kb/s)
* avg_CTCF_attach_time (averaged CTCF residence time, in s)
* proba_CTCF_occupancy (probability that a CTCF bound when a cohesin reaches a binding site)
* deltaT (time step)
* totalAreaMb (size of the simulated polymer)
The following methods allows then to simulate loop extrusion
* addCTCFsite() allows to add a CTCF site at a position with a given probability to block cohesin
* setAnchors() allows to set the position of two TAD anchors. The simulation will compute the distance between these anchors at each step (deltaT)
* run() allows to run the simulation for a given number of steps
* computeAnchor2AnchorDistance() allows to compute the distance between the two anchors
* computeAndGet_HiC_map() allows to compute the HiC map from these 1D simulations
