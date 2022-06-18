=============
Release notes
=============


Development version
===================

* Extended the Slater-Koster tabulation procedure to *f* electrons.
  If one of the two atoms contains *f* electrons in the valence,
  the SKF files are now written in the 'extended' format instead
  of the 'simple' format.

* Maintenance related to changes in ASE, Matplotlib and LibXC.


Version 1.0
===========

26 December 2019

* The GPAW-based atomic DFT calculator has been removed, as it was
  simply no longer needed (the native calculator now runs as
  fast and is easier to extend with more functionality).

* The native atomic DFT calculator has been renamed to 'AtomicDFT'
  and now resides in hotcent.atomic_dft.py. Its 'xcname' keyword
  argument (for the exchange-correlation functional) has furthermore
  been shortened to 'xc'.

* Calculation of Hubbard values (see examples/hubbard.py).

* Optional C-extensions allow for significantly faster atomic
  atomic DFT calculations (8-10x speedup) and construction of
  Slater-Koster tables (2-3x speedup).

* A new `stride` option for SlaterKosterTable.run() with default = 1.
  Setting higher integer values for `stride` means that the
  H and S integrals are only be explicitly calculated every
  `stride` points and will then be mapped on the final, denser grid
  using cubic spline interpolation. Use with care.

* No more backwards compatibility with Python2 (only Python3).

* The SlaterKosterTable.run() method now requires rmin, dr and N
  as arguments to specify the interatomic distances for which the
  Slater-Koster integrals are tabulated (rmin being the minimal
  distance, dr the grid spacing dr, and N the number of grid points).

* Major restructuring of hotcent.tools for fitting the confinement
  parameters, with also a new tutorial on this topic (`tutorial_2.ipynb`).


Version 0.9
===========

03 December 2019

* Start of versioning.
