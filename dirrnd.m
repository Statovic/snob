% Generate Dirichlet random numbers 
function x = dirrnd(params, n)

% shape = params; scale = 1
y = gamrnd(repmat(params(:)',n,1),1);
x = bsxfun(@rdivide, y, sum(y,2));

end