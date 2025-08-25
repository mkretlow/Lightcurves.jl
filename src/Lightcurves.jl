"""
### Purpose:
ALCTools.jl is a package for handling and analysing Asteroid Light Curves and Phase Curves.

### Author(s):
Mike Kretlow <mike@kretlow.de> [github.com/mkretlow]
"""
module Lightcurves

using ADTK
using AstroFun
using CSV
using DataFrames
using Interpolations
using LombScargle
using LinearAlgebra: norm
using LsqFit
#using Measures
#using PDM
using Plots
using Printf
using ProgressBars
using Statistics
using SQLite
using VSOP87

const version = v"0.1.0"

include("lightcurves.jl")
include("phasecurves.jl")

end # module
