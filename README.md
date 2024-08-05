# cosmicBallet

**cosmicBallet** is a python based N-Body Problem solver that has the capabilities to solve the N-Body Problem using various solvers that can be user specified depending on the required accuracy and runtime complexities.

## Installation Guide

This guide provides a simple installation of the package using specification files to create a clone of an enivronment as setup during development (i.e., no clashes or instabilities).

### Clone the Repository

Clone the main branch of the repository using the following command:

```
$ git clone --single-branch --branch main --depth 1 https://github.com/Sats2/cosmicBallet.git
```

And go to the directory where the repository clone exists (with cd command for instance).

### Package Installation

Create a virtual environment using [Anaconda](https://docs.anaconda.com/) (or [Miniconda](https://docs.anaconda.com/miniconda/) / [Micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for minimal installers). If unfimiliar with environment managers, refer to the help section linked here. 

1. For users on Windows/MacOS, simply copy and paste the following command in a conda (or mamba) terminal:
   For users on Windows or MacOS, the Binary System Merger simulations have not been tested. The other functionalities work as intended. Simply install the package within conda using -
   ```
   $ conda create --name cosmicBallet --file spec-file.txt
   $ python3 -m pip install .
   ```
3. For user on Unix, a full installation with Binary Merger Systems can be installed with either:
   
   a. Manual Installation - Use the following commands for installation
   ```
   $ conda create --name cosmicBallet --file spec-file-optional.txt
   $ pip install NRSur7dq2
   $ python3 -m pip install .
   ```
   b. Automatic Installation - Simply run the installation shell script (ensure the shell script has necessary permissions)
   ```
   $ source ./install.sh
   ```
   Permissions to the install script can be given with (chmod +x install.sh or chmod +755 install.sh)
   
   c. For minimal installation (without Binary Merger Systems):
   ```
   $ conda create --name cosmicBallet --file spec-file.txt
   $ python3 -m pip install .
   ```
The successful installation of the optional libraries are not guaranteed due to system requirements and has only been tested on WSL with Ubuntu 22.04LTS.

### Install Verification:

A successful installation can be verified using:

```
$ pip show cosmicBallet
```
which should show an output as follows:
```
Name: cosmicBallet
Version: 1.0.0
Summary: A Python based N-Body Solver
Home-page: https://github.com/Sats2/cosmicBallet.git
Author: Sathyamurthy Hegde
Author-email: sathyamhegde@outlook.com
License:
Location: /home/sathya/miniconda3/envs/cosmicball/lib/python3.11/site-packages
Requires: peppercorn
Required-by:
```

## Alternate Installation

The previously described installation method is preferred to prevent any dependency clashes causing unstable enstable environment conditions as not all libraries are available with conda installation (environment stability is not guaranteed with this method).

### Clone the Repository

Clone the main branch of the repository using the following command:

```
$ git clone --single-branch --branch main --depth 1 https://github.com/Sats2/cosmicBallet.git
```

### Create an Environment

Create a virtual environment with conda using:
```
$ conda create -f environment.yml
```

### Activate Environment

Activate the enivronment upon successful creation

```
$ conda activate cosmicBallet
```

### Install the Libraries

All essential libraries can be simply installed with pip. If the OS is Linux (or Windows Sub-system for Linux (WSL)) and a binary system merger visualization is necessary, the optional requirements need to be installed as well. The installation commands are:

```
$ pip install -r requirements.txt
```

and for the optional libraries
```
$ pip install -r requirements_optional.txt
```

### Package Installation

After the library installation, the package can be installed with:

```
$ python3 -m pip install .
```

## Usage

Refer to the tests folder of the repository for implementations.

## Documentation

The wiki (to do) within GitHub provides a useful guide to the package. Click here for the wiki.

## License

**_cosmicBallet_** is a python package created by Sathyamurthy Hegde (GitHub User account: Sats2) for the N-Body Problems with an MIT License.
