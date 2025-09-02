# MUSEgeom
A toolkit for geometry modeling for environmental domains in geoscience.

## Clone
We provide the commands to install MUSE. 
The repository includes some submodules necessary to make the code work. Please, clone it recursively:

- Clone recursively the repository into your local machine:
```
git clone --recursive git@github.com:mariannamiola/MUSEgeom.git
```
- If some submodules are not clone/update, please use git command in the root directory ROOT (where this README lies):
```
cd ${ROOT}
git submodule update --init --recursive
```

## Content of the repository
- `src`: source code
- `include`: functionalities library to support source codes
- `external`: external libraries
- `examples`: application data (input data and output metadata)
- `CMakeLists.txt`: build configuration file

## Dependences
MUSE code has some mandatory dependences (included as a submodule in _${ROOT}/external_):

- to manage command line arguments and options: [tclap](https://tclap.sourceforge.net/) (to clone and include in _${ROOT}/external_);
- to metadata the computational process: `cereal` (included as a submodule in _${ROOT}/external_).
- to geostatistics and stochastic computation (geostatslib);
- to manage geospatial data (GDAL, PROJ);
- to process polygonal/polyhedral meshes and offers geometric processing tools (cinolib, libigl, Triangle, Tetgen, fTetWild).

## Building
To build EWoPe source code, use the following pipeline:

```
cd ${ROOT}
mkdir build
cd build
cmake --DCMAKE_BUILD_TYPE=Release ..
cmake --build . --config Release
```

Executables will be made available in _${ROOT}/bin_ folder.

## Documentation
The complete set of options for each MUSE module is available [here](MUSE_Manual.pdf). 
MUSE can be easily built (see __Building__ section).


## Examples
To guarantee replicability, data of our examples are provided in the _example_ folder. 


## Contributors
- Marianna Miola (CNR-IMATI, Genova, Italy), email: marianna.miola@cnr.it

## Citing us
If you use our tool in your academic projects, please consider citing us using the following BibTeX entry:
```

```

## Acknowledgment