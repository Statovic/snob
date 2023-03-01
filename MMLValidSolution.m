% Find valid root with minimum codelength 
function bestX = MMLValidSolution(N, K, J, delta, x0)

L = length(x0);
minMsglen = +inf; bestX = delta(J);
for i = 1:L
    m = pca_codelength(N, K, J, x0(i), delta);
    if(isreal(m) && isreal(x0(i)))
        if(m < minMsglen)
            minMsglen = m; 
            bestX = x0(i);
        end
    end
end

end
