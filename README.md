# MA473-Computational-Finance
Lab Assignments of MA473 - Computational Finance 2024.

# Some missing things
1. Lab 01 - The 2nd question will be solved by discretising the pde and applying implicit and explicit methods as learnt before in scientific computing.

2. Lab 03 - The Theta greek can be calculated just like gamma and delta in the plots.m file. However for vega and rho, the schemes will be called multiple time with the variation in the required parameter, then the difference of matrix obtained divided by change in parameter will give the value of those greeks over different time and space points.

3. Lab 05 - The plots for order of convergence and delx vs max Error plots were giving absurd results when calculated by giving finite limit to the integral, and NaN and INF values if limit was passed to be infinite.

<TO BE UPDATED>