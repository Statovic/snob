% Generate random numbers from the Pareto Type IV distribution
function x = prtrnd(n, sigma, gam, alpha)

u = rand(n,1);
x = sigma * ((1 - u).^(-1/alpha) - 1).^(gam);

end