[![Build Status](https://travis-ci.org/ibex-team/ibex-lib.svg?branch=master)](https://travis-ci.org/ibex-team/ibex-lib)
[![Build status](https://ci.appveyor.com/api/projects/status/9w1wxhvymsohs4gr/branch/master?svg=true)](https://ci.appveyor.com/project/Jordan08/ibex-lib-q0c47/branch/master)

Linearization methods: AbsTaylor & TrustRegions
=======

<h1>AbsTaylor</h1>
<p>This method is based on a first order Taylor-based linearization technique, which can be used for finding feasible solutions of Numerical Constraint Optimization Problems (NCOPs).</p>
<p>Unlike the previous proposal, X-Taylor, the expansion point of the Taylor form is not limited to the bounds of the domains, thus producing, in general, a tighter approximation.</p>
<p>In order to apply the strategy the following steps must be follow:</p>
<ul>
<li>
<p>Decompose the interval partial derivatives into two values: <img src="https://i.upmath.me/svg/J_%7Bij%7D%3Dc_%7Bij%7D%2B%5Cmathbf%7Bd_%7Bij%7D%7D" alt="J_{ij}=c_{ij}+\mathbf{d_{ij}}" />, where the value <img src="https://i.upmath.me/svg/c_%7Bij%7D" alt="c_{ij}" /> is the midpoint of the interval derivative.</p>
</li>
<li>
<p>Then, we propose the following Taylor-based inner linearization:
<img src="https://i.upmath.me/svg/%20g_j(x)%20%5Cleq%20h_j(x)%3Dg_j(x')%20%2B%20%5Csum_i%5En%20%5Cleft(c_%7Bij%7D%20%5Ccdot%20(x_i%20-%20x_i')%20%2B%20%7Cd_%7Bij%7D%7C%0A%20%20%5Ccdot%20%7Cx_i%20-%20x_i'%7C%20%5Cright)" alt=" g_j(x) \leq h_j(x)=g_j(x') + \sum_i^n \left(c_{ij} \cdot (x_i - x_i') + |d_{ij}|
  \cdot |x_i - x_i'| \right)" /></p>
</li>
</ul>
<p>where <img src="https://i.upmath.me/svg/x'%5Cin%5Cmathbb%7BR%7D%5En" alt="x'\in\mathbb{R}^n" />
can be any point inside the box x, and <img src="https://i.upmath.me/svg/%7Cd_%7Bij%7D%7C%3D%5Coverline%7Bd_%7Bij%7D%7D%3D%5Cunderline%7Bd_%7Bij%7D%7D" alt="|d_{ij}|=\overline{d_{ij}}=\underline{d_{ij}}" /></p>
<ul>
<li>For finding a promising feasible solution in the inner region, we propose to minimize a linearization of the objective function subject to the constraint system generated before:</li>
</ul>
<p align="center"><img align="center" src="https://i.upmath.me/svg/%0A%20%20%20%20%5Cbegin%7Btabular%7D%7Bl%20l%7D%0A%20%20%20%20%20%20%20%20%24%5Cmin%5Climits_%7Bx%20%5Cin%20%7B%5Cbf%20x%7D%7D%24%20%26%20%24f(x')%20%2B%20%5Csum%5Climits_i%5En%20%5Cmbox%7Bmid%7D(%7B%5Cbf%20J_%7Bij%7D%7D)%0A%20%20%20%20%20%20%20%20%5Ccdot%20(x_i%20-%20x_i')%24%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20s.t%20%26%20%24g_j(x')%20%2B%20%5Csum%5Climits_i%5En%20%5Cleft(c_%7Bij%7D%20%5Ccdot%20(x_i%20-%20x_i')%20%2B%20%7Cd_%7Bij%7D%7C%0A%20%20%20%20%20%20%20%20%20%20%20%20%5Ccdot%20%7Cx_i%20-%20x_i'%7C%20%5Cright)%20%5Cleq%200%24.%20%5C%5C%0A%20%20%20%20%5Cend%7Btabular%7D%0A" alt="
    \begin{tabular}{l l}
        $\min\limits_{x \in {\bf x}}$ &amp; $f(x') + \sum\limits_i^n \mbox{mid}({\bf J_{ij}})
        \cdot (x_i - x_i')$ \\
            s.t &amp; $g_j(x') + \sum\limits_i^n \left(c_{ij} \cdot (x_i - x_i') + |d_{ij}|
            \cdot |x_i - x_i'| \right) \leq 0$. \\
    \end{tabular}
" /></p>
<ul>
<li>In order to deal with absolute values functions, we replace the expresions <img src="https://i.upmath.me/svg/%7Cx_%7Bi%7D-x_%7Bi%7D'%7C" alt="|x_{i}-x_{i}'|" /> by auxiliary variables <img src="https://i.upmath.me/svg/u_i" alt="u_i" />. We also add the constraints <img src="https://i.upmath.me/svg/u_i%20%5Cgeq%20x_i%20-%20x_i'%24%20and%20%24u_i%20%5Cgeq%20-(x_i%20-%20x_i')" alt="u_i \geq x_i - x_i'$ and $u_i \geq -(x_i - x_i')" />, because the partial derivative of the constraint functions is positive. Finally we solve the new system using the Simplex method:</li>
</ul>
<p align="center"><img align="center" src="https://i.upmath.me/svg/%0A%5Cbegin%7Btabular%7D%7Bl%20l%7D%0A%20%20%20%20%20%20%20%20%24%5Cmin%5Climits_%7Bx%20%5Cin%20%7B%5Cbf%20x%7D%7D%24%20%26%20%24f(x')%20%2B%20%5Csum%5Climits_i%5En%20%5Cmbox%7Bmid%7D(%7B%5Cbf%20J_%7Bij%7D%7D)%20%0A%20%20%20%20%20%20%20%20%5Ccdot%20(x_i%20-%20x_i')%24%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20s.t%20%26%20%24g_j(x')%20%2B%20%5Csum%5Climits_i%5En%20%5Cleft(c_%7Bij%7D%20%5Ccdot%20(x_i%20-%20x_i')%20%2B%20%7Cd_%7Bij%7D%7C%20%5Ccdot%20u_i%20%5Cright)%20%5Cleq%200%24%2C%20~~%24%5Cforall%20j%3D1%5Cldots%20m%24%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%20%24u_i%20%5Cgeq%20x_i%20-%20x_i'%24%2C%20~~%24%5Cforall%20i%3D1%5Cldots%20n%24%20%5C%5C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%26%20%20%24u_i%20%5Cgeq%20-(x_i%20-%20x_i')%24%2C%20~~%24%5Cforall%20i%3D1%5Cldots%20n%24.%20%0A%20%20%20%20%5Cend%7Btabular%7D%0A" alt="
\begin{tabular}{l l}
        $\min\limits_{x \in {\bf x}}$ &amp; $f(x') + \sum\limits_i^n \mbox{mid}({\bf J_{ij}}) 
        \cdot (x_i - x_i')$ \\
            s.t &amp; $g_j(x') + \sum\limits_i^n \left(c_{ij} \cdot (x_i - x_i') + |d_{ij}| \cdot u_i \right) \leq 0$, ~~$\forall j=1\ldots m$ \\
             &amp;  $u_i \geq x_i - x_i'$, ~~$\forall i=1\ldots n$ \\
             &amp;  $u_i \geq -(x_i - x_i')$, ~~$\forall i=1\ldots n$. 
    \end{tabular}
" /></p>
<h2>Example</h2>
<p>Let consider the following system of constraints:</p>
<p align="center"><img align="center" src="https://i.upmath.me/svg/%0A%5Cbegin%7Beqnarray*%7D%0A%20%20%20%20g_%7B1%7D(x)%26%3D%26x_%7B1%7D%5E5%2B0.5%5Ccos(x_%7B1%7D)%2B%5Csin(x_%7B2%7D)-2x_%7B2%7D-0.2%20%5Cleq%200%20%5C%5C%0A%20%20%20%20g_%7B2%7D(x)%26%3D%26-x_%7B1%7D%2Bx_%7B2%7D%5E%7B2%7D-1%20%5Cleq%200%20%5C%5C%0A%5Cend%7Beqnarray*%7D%0A" alt="
\begin{eqnarray*}
    g_{1}(x)&amp;=&amp;x_{1}^5+0.5\cos(x_{1})+\sin(x_{2})-2x_{2}-0.2 \leq 0 \\
    g_{2}(x)&amp;=&amp;-x_{1}+x_{2}^{2}-1 \leq 0 \\
\end{eqnarray*}
" /></p>
<p>in the box <img src="https://i.upmath.me/svg/%7B%5Cbf%20x%7D%3D%5B-1%2C1%5D%5Ctimes%20%5B0%2C1%5D" alt="{\bf x}=[-1,1]\times [0,1]" />. By applying AbsTaylor, and using the midpoint of the box as the expansion point, we obtain the following linearization:</p>
<p align="center"><img align="center" src="https://i.upmath.me/svg/%0A%5Cbegin%7Beqnarray*%7D%0A%20%20%20%20g_%7B1%7D%5Cleq%26%202.5x_%7B1%7D%2B2.9207%7Cx_%7B1%7D%7C-1.2298x_%7B2%7D%2B0.2298%7Cx_%7B2%7D-0.5%7C%2B0.4944%20%26%5Cleq%200%5C%5C%0A%20%20%20%20g_%7B2%7D%5Cleq%26%20-x_%7B1%7D%2Bx_%7B2%7D%2B%7Cx_%7B2%7D-0.5%7C-1.25%20%26%5Cleq%200%5C%5C%0A%5Cend%7Beqnarray*%7D%0A" alt="
\begin{eqnarray*}
    g_{1}\leq&amp; 2.5x_{1}+2.9207|x_{1}|-1.2298x_{2}+0.2298|x_{2}-0.5|+0.4944 &amp;\leq 0\\
    g_{2}\leq&amp; -x_{1}+x_{2}+|x_{2}-0.5|-1.25 &amp;\leq 0\\
\end{eqnarray*}
" /></p>
<h1>IterativeTaylor</h1>
<p>IterativeTaylor is a simple iterative method for finding promising feasible solutions in a box <img src="https://i.upmath.me/svg/%7B%5Cbf%20x%7D" alt="{\bf x}" />.</p>
<p>Starting from an initial point <img src="https://i.upmath.me/svg/x_p" alt="x_p" /> (e.g., the midpoint
of the box <img src="https://i.upmath.me/svg/%7B%5Cbf%20x%7D" alt="{\bf x}" />), IterativeTaylor first attempts to find a candidate feasible solution <img src="https://i.upmath.me/svg/x_c" alt="x_c" /> in <img src="https://i.upmath.me/svg/%7B%5Cbf%20x%7D" alt="{\bf x}" /> by using
any LoupFinder (such as AbsTaylor or XTaylor) with <img src="https://i.upmath.me/svg/x_c" alt="x_c" /> as the expansion point of the linearization.</p>
<p>If it is successful, then the algorithm performs a loop while the candidate solution is improved. In each iteration, the size of <img src="https://i.upmath.me/svg/%7B%5Cbf%20x%7D" alt="{\bf x}" /> is reduced and a new feasible solution <img src="https://i.upmath.me/svg/x_c" alt="x_c" /> is searched in the new box by applying a linearization technique (or a set of linearization) and the simplex method.</p>
<h3>Authors: Victor Reyes and Ignacio Araya</h3>
<h5>For citing use the paper: Reyes, V., &amp; Araya, I. (2020). AbsTaylor: upper bounding with inner regions in nonlinear continuous global optimization problems. Journal of Global Optimization, 1-17.</h5>

