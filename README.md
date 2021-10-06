# NoRCB
Non-orthogonal recursive coordinate bisection is a generalization of recursive coordinate bisection for lagrangian problem that require partitioning of particles among processing elements. It is essentially a C++ class that performs several operations on a set of elements that are distributed in space.
This library is meant to be added as a git submodule to other projects.

# Dependencies
- CGAL
- MPI

# Build 
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j
```

# How to
