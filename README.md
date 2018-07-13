# aN2MD v0.1

aN2MD is a program designed to compare NMR data with simulations and structures. 

It has two main utilities:
* Creation of .itp files availables for GROMACS from two .str files. 
* NOE distances violation analysis of GROMACS simulations or single structures.



## Getting Started

These instructions will get you a copy of the project up and running on your
local machine for development and testing purposes.

## Prerequisites 

This program is developped with Python3 and designed for linux (Ubuntu 16.04+).
This program requires the following Python modules : 
* [mdtraj](http://mdtraj.org/1.8.0/installation.html)
* [numpy](https://scipy.org/install.html)
* [pandas](https://pandas.pydata.org/)
* [matplotlib](https://matplotlib.org/users/installing.html)


## Installing

Here is a step by step guide to install and run the program. Before starting,
check the Prerequesites section and ensure that all the required modules are
correctly installed. 

1. Clone the github directory (requires git).

```
mkdir aN2MD
git clone https://github.com/PaulaMilanRguez/aN2MD.git aN2MD
```

No further installation is required to run the program. 


## Running the Demo

Run the program Main_NCp7.py to test the installation 

```
python3 Main_NCp7.py
```


## Changelog

Version 0.1 :
* First deployable version, uploaded on github 

## Authors

* Paula Milan Rodriguez
* Marco Pasi

## License

Program under Gnu General Public License v2.0

## Acknowledgments 

