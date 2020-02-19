# CompactBases.jl

A package for representing various bases constructed from basis
functions with compact support as quasi-arrays.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaApproximation.github.io/CompactBases.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaApproximation.github.io/CompactBases.jl/dev)
[![Build Status](https://travis-ci.com/JuliaApproximation/CompactBases.jl.svg?branch=master)](https://travis-ci.com/JuliaApproximation/CompactBases.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/JuliaApproximation/CompactBases.jl?svg=true)](https://ci.appveyor.com/project/JuliaApproximation/CompactBases-jl)
[![Codecov](https://codecov.io/gh/JuliaApproximation/CompactBases.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaApproximation/CompactBases.jl)

This package implements bases with compactly supported bases functions
as quasi-arrays where one one axes is continuous and the other axis is
discrete (countably infinite), as implemented in
[QuasiArrays.jl](https://github.com/JuliaApproximation/QuasiArrays.jl)
and
[ContinuumArrays.jl](https://github.com/JuliaApproximation/ContinuumArrays.jl). The
bases available initially are various finite-differences (FD),
finite-elements discrete-variable representation (FE-DVR), and
B-splines, with the possibility of adding more in the future.
