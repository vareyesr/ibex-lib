[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

Linearization methods: AbsTaylor & TrustRegions
=======

# AbsTaylor

This method is based on a first order Taylor-based linearization technique, which can be used for finding feasible solutions of Numerical Constraint Optimization Problems (NCOPs).

Unlike the previous proposal, X-Taylor, the expansion point of the Taylor form is not limited to the bounds of the domains, thus producing, in general, a tighter approximation.

In order to apply the strategy the following steps must be follow:

* Decompose the interval partial derivatives into two values: <a href="https://www.codecogs.com/eqnedit.php?latex=J_{ij}=c_{ij}&plus;\mathbf{d_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?J_{ij}=c_{ij}&plus;\mathbf{d_{ij}}" title="J_{ij}=c_{ij}+\mathbf{d_{ij}}" /></a>, where the value c_{ij} is the midpoint of the interval derivative.

* Then, we propose the following Taylor-based inner linearization:
<center>
<a href="https://www.codecogs.com/eqnedit.php?latex=g_{j}\leq&space;h_{j}&space;=&space;g_{j}(x')&plus;\sum_{i}^{n}(c_{ij}\cdot&space;(x_{i}-x_{i}')&plus;|d_{ij}|\cdot&space;|x_{i}-x_{i}'|&space;)" target="_blank"><img src="https://latex.codecogs.com/svg.latex?g_{j}\leq&space;h_{j}&space;=&space;g_{j}(x')&plus;\sum_{i}^{n}(c_{ij}\cdot&space;(x_{i}-x_{i}')&plus;|d_{ij}|\cdot&space;|x_{i}-x_{i}'|&space;)" title="g_{j}\leq h_{j} = g_{j}(x')+\sum_{i}^{n}(c_{ij}\cdot (x_{i}-x_{i}')+|d_{ij}|\cdot |x_{i}-x_{i}'| )" /></a>	
</center>
where <a href="https://www.codecogs.com/eqnedit.php?latex=x'\in&space;\mathbb{R}^n" target="_blank"><img src="https://latex.codecogs.com/svg.latex?x'\in&space;\mathbb{R}^n" title="x'\in \mathbb{R}^n" /></a> 
can be any point inside the box x, and <a href="https://www.codecogs.com/eqnedit.php?latex=|d_{ij}|=\overline{d_{ij}}=\underline{d_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?|d_{ij}|=\overline{d_{ij}}=\underline{d_{ij}}" title="|d_{ij}|=\overline{d_{ij}}=\underline{d_{ij}}" /></a>.

* For finding a promising feasible solution in the inner region, we propose to minimize a linearization of the objective function subject to the constraint system generated before:
<center>
<a href="https://www.codecogs.com/eqnedit.php?latex=\min\limits_{x&space;\in&space;x}&space;f(x')&space;&plus;&space;\sum_i^n&space;mid(J_{ij})\cdot&space;(x_i&space;-&space;x_i')" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\min\limits_{x&space;\in&space;x}&space;f(x')&space;&plus;&space;\sum_i^n&space;mid(J_{ij})\cdot&space;(x_i&space;-&space;x_i')" title="\min\limits_{x \in x} f(x') + \sum_i^n mid(J_{ij})\cdot (x_i - x_i')" /></a> <br/>
</center>
<center>
<a href="https://www.codecogs.com/eqnedit.php?latex=s.t&space;\mbox{&space;}&space;g_j(x')&space;&plus;&space;\sum_i^n&space;\left(c_{ij}&space;\cdot&space;(x_i&space;-&space;x_i')&space;&plus;&space;|d_{ij}|\cdot&space;|x_i&space;-&space;x_i'|&space;\right)&space;\leq&space;0." target="_blank"><img src="https://latex.codecogs.com/svg.latex?s.t&space;\mbox{&space;}&space;g_j(x')&space;&plus;&space;\sum_i^n&space;\left(c_{ij}&space;\cdot&space;(x_i&space;-&space;x_i')&space;&plus;&space;|d_{ij}|\cdot&space;|x_i&space;-&space;x_i'|&space;\right)&space;\leq&space;0." title="s.t \mbox{ } g_j(x') + \sum_i^n \left(c_{ij} \cdot (x_i - x_i') + |d_{ij}|\cdot |x_i - x_i'| \right) \leq 0." /></a>
</center>

In order to deal with absolute values functions, we replace the expresions <a href="https://www.codecogs.com/eqnedit.php?latex=|x_{i}-x_{i}'|" target="_blank"><img src="https://latex.codecogs.com/svg.latex?|x_{i}-x_{i}'|" title="|x_{i}-x_{i}'|" /></a> by auxiliary variables <a href="https://www.codecogs.com/eqnedit.php?latex=u_i" target="_blank"><img src="https://latex.codecogs.com/svg.latex?u_i" title="u_i" /></a>. We also add the constraints <a href="https://www.codecogs.com/eqnedit.php?latex=u_i&space;\geq&space;x_i&space;-&space;x_i'" target="_blank"><img src="https://latex.codecogs.com/svg.latex?u_i&space;\geq&space;x_i&space;-&space;x_i'" title="u_i \geq x_i - x_i'" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=u_i&space;\geq&space;-(x_i&space;-&space;x_i')" target="_blank"><img src="https://latex.codecogs.com/svg.latex?u_i&space;\geq&space;-(x_i&space;-&space;x_i')" title="u_i \geq -(x_i - x_i')" /></a>, because the partial derivative of the constraint functions is positive. Finally we solve the new system using the Simplex method:

# TrustRegions








## Authors: Victor Reyes and Ignacio Araya

#### For citing use the paper: Reyes, V., & Araya, I. (2020). AbsTaylor: upper bounding with inner regions in nonlinear continuous global optimization problems. Journal of Global Optimization, 1-17.

