# cosmicBallet

**cosmicBallet** is a python based N-Body Problem solver that has the capabilities to solve the N-Body Problem using various solvers that can be user specified depending on the required accuracy and runtime complexities.

## Installation

### Clone the Repository

Clone the main branch of the repository using the following command:

```
git clone --single-branch --branch main --depth 1 https://github.com/Sats2/cosmicBallet.git
```

### Create an Environment

Create a virtual environment using [Anaconda](https://docs.anaconda.com/) (or [Miniconda](https://docs.anaconda.com/miniconda/) / [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for minimal installers). If unfimiliar with environment managers, refer to the help section linked here. 
```
conda create -f environment.yml
```

### Activate Environment

Activate the enivronment upon successful creation

```
conda activate cosmicBallet
```

### Install the Libraries

All essential libraries can be simply installed with pip. If the OS is Linux (or Windows Sub-system for Linux (WSL)) and a binary system merger visualization is necessary, the optional requirements need to be installed as well. The installation commands are:

```
pip install -r requirements.txt
```

and for the optional libraries
```
pip install -r requirements_optional.txt
```
The successful installation of the optional libraries are not guaranteed due to system requirements and has only been tested on WSL with Ubuntu 22.04LTS.

### Package Installation

After the library installation, the package can be installed with:

```
pip install .
```

## Usage

Refer to the tests folder of the repository for implementations.

## Documentation

The wiki within GitHub provides a useful guide to the package. Click here for the wiki.

## License

**_cosmicBallet_** is a python package created by Sathyamurthy Hegde (GitHub User account: Sats2) for the N-Body Problems with an MIT License.
