# A list of global TODO items

This list should contain items that are not (yet) linked to a particular part of
the code.

\todo **Important** Convert amrvacdef.t to a module

\todo **Important** Convert the physics modules to actual modules

\todo Allow runs with more CPUs than coarse grid blocks

\todo Implement the hlld solver

\todo Individual time step scheme

\todo Make tracers work with roe solvers

\todo Include the python doxygen documentation

\todo Include the IDL documentation

\todo Replace magic numbers with symbolic constants for clarity:

* neighbor_type
* errorestimate

\todo Replace strings with symbolic constants, which increases performance:

* typelimiter
* typeaxial
* typeB
* typefull1
* typeprolonglimit
* typeghostfill
* etc.

See [this benchmark](https://github.com/jannisteunissen/fortran_benchmarks) for
the performance difference.
