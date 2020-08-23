# Reversible-jump Change Points in Julia

This is reimplementation of an example given by Green in his paper on
reversible-jump Markov chain Monte Carlo, concerned with change-point detection.
The original implementation was done in Fortran and used Gibbs sampling for
within-model moves; this version is written in Julia and uses Hamiltonian Monte
Carlo instead.

Note: this was written before Julia 1.0 was released, and I haven't updated it
to work with that or newer releases. While the code here might be useful for
reference (or perhaps even use—who knows), the main files of interest are the
two notes in the `doc` directory: the [first](doc/rjcp.pdf) summarising the
approach, implementation, and results; and the [second](doc/green.pdf) something
of a tutorial/exposition on the theory behind reversible-jump Markov chain Monte
Carlo.