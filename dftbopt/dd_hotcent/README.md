# Hotcent

Calculating one- and two-center Slater-Koster integrals,
based on parts of the [Hotbit](https://github.com/pekkosk/hotbit/)
code. The development of Hotcent was started as part of the
following study:

M. Van den Bossche, J. Chem. Phys. A. **2019**, 123 (13), 3038-3045
[(doi)](https://dx.doi.org/10.1021/acs.jpca.9b00927).


## Features

* The main use for Hotcent is generating Slater-Koster tables,
typically in the [".skf"](https://www.dftb.org/fileadmin/DFTB/public/misc/slakoformat.pdf)
format to be used with DFTB codes such as DFTB+.

* The code allows to use different confinement potentials for
the different valence wave functions and for the electron density
(which determines the effective potential used in the Hamiltonian
integrals).

* Both the potential superposition and density superposition
schemes are available.

* With regards to exchange-correlation functionals, the PW92
(LDA) functional is natively available, and other LDA/GGA
functionals can be applied through integration with the PyLibXC
module shipped with [LibXC](https://www.tddft.org/programs/libxc).
Hybrid and meta-GGA functionals cannot currently be used in
Hotcent.


## Installation

* Set up the [ASE](https://wiki.fysik.dtu.dk/ase/) Python module
  (version 3.21.1 or newer).

* Clone / download the Hotcent repository and update the
`$PYTHONPATH` accordingly, e.g. like this:
```shell
cd <your_location_of_choice>
git clone https://gitlab.com/mvdb/hotcent
export PYTHONPATH=$PWD/hotcent:$PYTHONPATH
```

* If you want significantly faster calculations (who doesn't?),
you will want to build the optional C-extensions:
```shell
python setup.py build_ext --inplace
```
If you have Cython installed and want/need to regenerate the C-code,
just run:
```shell
python setup.py build_ext --inplace --use-cython
```

* If you want to use functionals other than the PW92 LDA (again, who doesn't?),
the [PyLibXC](https://www.tddft.org/programs/libxc/installation/#python-library)
module needs to be available, which provides a Python interface to all
LibXC functionals. A recent LibXC version is required (>= 4.3.4).
Installing this module can e.g. be done as follows (modify as needed):
```shell
cd <your_libxc_directory>
export PYTHONPATH=$PWD/lib/python3.8/site-packages/:$PYTHONPATH
python setup.py install --prefix=$PWD
```
