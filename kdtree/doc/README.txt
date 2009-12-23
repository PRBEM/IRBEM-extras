kdtree - a C library for building and accessing k-dimensional
trees. These are handy for finding nearest neighbors in a Pythagorean
or maximum coordinate distance sense.

**** Installation ****
1. See the Makefile on how to compile a .dll, .so, or .oct
file. Note: this implementation supports OpenMP parallelization

2. For Matlab/Octave users, you'll need to add two folders to your
path: .../kdtree/matlab first, then .../kdtree. If you accidentally
put the matlab folder second, Matlab/Octave will try to load the
.dll/.so file as a MEX file when you call the "kdtree" function, and
this won't work. In a later rev, I'll rename the .dll and .so
something else (like libkdtree.*) so it won't have this name conflict.

**** Algorithm ****

The kdtree search is approximately log(n) for finding k nearest
neighbors in a point cloud of n points. More at
http://en.wikipedia.org/wiki/Kd-tree

**** Files of interest ****
Makefile - has make commands that work on CentOS, cygwin, and MinGW

kdree.h - C function prototypes

matlab/kdtree.m - Matlab/Octave wrapper

matlab/kdtree_demo.m - Matlab/Octave demo

Paul O'Brien
paul.obrien@aero.org
