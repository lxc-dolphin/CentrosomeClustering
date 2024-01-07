# Stochastic Simulation for [centrosome movement and spindle assembly](https://www.molbiolcell.org/doi/full/10.1091/mbc.E22-10-0485)
## Features
- Using Langevin equation to simulate centrosome dynamics and stochastic fluctuations.


    $\ \frac{d\bold{X}}{dt} = 1/\gamma (\bold{F}_{cs}+\bold{F}_{rad}) + \eta(t) \$


- Integrating inter-CS force energy and radial energy potentials.
- Using step-adaptive method to solve stochastic differential equations. 
  
## Requirements
- [Julia](https://julialang.org/)
- Julia Packages:
  LinearAlgebra, DifferentialEquations, Distances, Combinatorics, Plots, Clustering, Interpolations.

## Example


<p align="center">
    <br>
    <hr />
    <img src="https://github.com/lxc-dolphin/CentrosomeClustering/blob/main/sup/CS_clustering1.gif" width="800">
    <hr />  
    <br>
<p>

