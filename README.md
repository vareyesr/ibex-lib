[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

Linearization methods: AbsTaylor & TrustRegions
=======

# AbsTaylor

This method is based on a first order Taylor-based linearization technique, which can be used for finding feasible solutions of Numerical Constraint Optimization Problems (NCOPs).

Unlike the previous proposal, X-Taylor, the expansion point of the Taylor form is not limited to the bounds of the domains, thus producing, in general, a tighter approximation.

In order to apply the strategy the following steps must be follow:

* Decompose the interval partial derivatives into two values: <a href="https://www.codecogs.com/eqnedit.php?latex=J_{ij}=c_{ij}&plus;\mathbf{d_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?J_{ij}=c_{ij}&plus;\mathbf{d_{ij}}" title="J_{ij}=c_{ij}+\mathbf{d_{ij}}" /></a>, where the value c_{ij} 


# TrustRegions








## Authors: Victor Reyes and Ignacio Araya

#### For citing use the paper: Reyes, V., & Araya, I. (2020). AbsTaylor: upper bounding with inner regions in nonlinear continuous global optimization problems. Journal of Global Optimization, 1-17.

