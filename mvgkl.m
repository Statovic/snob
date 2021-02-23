%MVGKL    Kullback-Leibler divergence between two multivariate Gaussians.
%    KL = MVGKL(m1,S1,m2,S2) returns the KL divergence between two
%    multivariate Gaussian distributions P1 and P2. P1 has parameters m1
%    (mean) and S1 (covariance matrix). P2 has parameters M2 (mean) and S2
%    (covariance matrix).
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function kl = mvgkl(m1, S1, m2, S2)

%% error checking
if(~iscolumn(m1) || ~iscolumn(m2))
    error('mean parameters must be column vectors');
end
if(length(m1) ~= length(m2))
    error('mean vectors are different dimenions');
end
if(~ismatrix(S1) || ~ismatrix(S2))
    error('covariance parameters must be matrices');
end
if(~all(size(S1) == size(S2)))
    error('covariance parameters are different dimenions');
end

% d-variate Gaussian
d = length(m1);

[R1,P1] = cholcov(S1,0); % Cholesky decomposition of covariance matrices
[R2,P2] = cholcov(S2,0);

if(any([P1,P2]))
    error('covariance matrices are not positive definite');
end

%% Compute KL divergence
sqTerm = sum( ((m2-m1)' / R2).^2 ); % Squared term

logDetS1 = 2*sum(log(diag(R1)));    % log |S1|
logDetS2 = 2*sum(log(diag(R2)));    % log |S2|

% KL divergence
kl = trace(R2\(R2'\S1)) + sqTerm - d + logDetS2 - logDetS1;
kl = kl / 2;


end
