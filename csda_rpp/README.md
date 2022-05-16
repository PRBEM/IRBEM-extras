# csda_rpp
SEE rate calculations using CSDA and RPP approximations
SEE - Single event effects
CSDA - Continuous slowing down approximation
RPP - Right parallel piped

This software is provided under the terms of the license provided in LICENSE.txt.

Initial coder: Paul O'Brien (paul.obrien@aero.org)
Contributors: Matt Worstell, Scott Davis

The following python files are included. Documentation is provided via docstrings. 
- degraded_spectra.py - compute CSDA degraded spectra for idealized shielding 
    geometries. Includes a demo.
- rpp.py - compute SEE rate in parts using RPP approximation. Includes a demo
    and path length distribution plots.
- shielded_part.py - wrapper with extensive caching for proton and ion (RPP) 
    part SEE calculations, including degraded spectrum computations
    wraps shielded_part.py and rpp.py. Includes demo of calling/caching.
- interpolating_matrix.py - compute matrix that interpolates between two 1-D 
    grids. Includes a demo
- creme.py - read CREME96 output files for ingest/comparisons. Includes a demo.
- demo.py - a demonstration of various package capabilities, comparisons to CREME

imported standard/common python modules:
numpy, matplotlib, os, sys, warnings, scipy, datetime, inspect

data included in repository:

Range-Energy Tables: (data/ sub-folder)
- elements.csv - table of atomic elements, Z's, A's, and MAI's
- NISTRangeData_protons.csv - NIST range table for protons
- ScreamRangeData_protons.csv - SCREAM range table for protons
- ScreamRangeData_electrons.csv  - SCREAM range table for electrons
- SRIM_Al.csv - SRIM range table in Al
- SRIM_Si.csv - SRIM range table in Si

Sample CREME files for comparison: (data/ sub-folder)
- q800a90i100.DLT - differential LET spectrum
- q800a90i100Z2.DLT - differential LET spectrum Z>=2
- q800a90i.FLX - incident ion energy spectrum
- q800a90i100.LET - integral LET spectrum
- q800a90i100.TFX - degraded ion energy spectrum
- q800a90i100.PUP - proton part SEE rates
- q800a90i100Z2.HUP - ion part SEE rates

Figures: (figures/ sub-folder)
- creme_fluxes.png - incident and degraded fluxes from CREME
- csda_fluxes.png - incident and degraded fluxes from CREME and CSDA
- let_spectrum.png - LET spectra from CREME and CSDA
- LB88_fig1.png - Chord length distribution for Figure 1 from Luke and Buehler 1988
- LB88_fig2.png - Chord length distribution for Figure 2 from Luke and Buehler 1988
- SCREAM_degraded_slow_fast_slab.png - comparison of SCREAM degraded fluxes for slabs
- SRIM_degraded_slow_fast.png - comparison of SRIM degraded fluxes
- SRIM_let_slow_fast.png - SRM LET spectrum
- NIST_degraded_slow_fast_SCREAM_slab.png - NIST/SCREAM degraded slab fluxes
- SCREAM_SRIM_degraded_slow_fast_sphere.png - SCREAM/SRIM degraded sphere fluxes
- SCREAM_electrons_degraded_slab.png - electron CSDA degraded flux

SCREAM data files provided by Scott Messenger (thank you, Scott!)
SRIM can be found at srim.org
CREME can be found at creme.isde.vanderbilt.edu

The calculations performed here can be found in ATR-2022-00960, 
Using Continuous Slowing Down Approximation and the Right Parallel 
Piped Model to Estimate Single Event Effects Rates by T.P. O'Brien
available from library.mailbox@aero.org
