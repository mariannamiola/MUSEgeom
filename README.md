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
To build TOOL source code, use the provided build script, which supports two options:

### Build all (GDAL + PROJ + TOOL)
This will compile the external dependencies (PROJ and GDAL), and then build the TOOL application:
```
cd ${ROOT}
./build.sh
```

### Build only the TOOL (if dependencies are already built)
Use this option if PROJ and GDAL have already been built and installed.
```
cd ${ROOT}
./build.sh --only-tool
```

Executables will be made available in _${ROOT}/bin_ folder.

## Documentation
The tool can be easily built (see __Building__ section).
To get an overview of all available command-line options:

```
./tool --help
```

To display help specific to each functionalities:
```
./tool --help-vol
```

## Examples
To guarantee replicability, data of our examples are provided in the _example_ folder. 


## Contributors
- Marianna Miola (CNR-IMATI, Genova, Italy), email: marianna.miola@cnr.it

## Citing us
If you use our tool in your academic projects, please consider citing us using the following BibTeX entry:
```

```

## Acknowledgment
