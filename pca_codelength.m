% MML codelength given an estimate of s2
function [msglen, assertion, detail] = pca_codelength(N, K, J, s2, delta)

alpha = sqrt(delta(1:J) - s2);

D = J*K - J*(J+1)/2;
params = 1 + J + D;

%% assertion
a2 = alpha .* alpha;
h_s = log(s2)/2;
h_R = J*log(2) + (K*J)/2*log(pi) - logmvgamma(J, K/2);
h_a = -J*log(2) - (J^2/2)*log(pi) - (J*J)/2*log(s2) + logmvgamma(J, J/2) + logmvbeta(J, K/2, J/2) ...
    -sum((K-J) * log(alpha)) + sum( (K+J)/2 * log(s2 + a2) );
h = h_s + h_R + h_a - gammaln(J+1);
F = params/2*log(N) + (J+1)/2*log(2) + 1/2*log(K-J) - (J*(K-J)+1)/2*log(s2) ...
    +1/2*sum((4*(K-J)+2)*log(alpha)) - (1/2)*sum((K+1)*log(a2 + s2));
assertion = h + F + (mml_const(params) - params/2);

%% detail
L = N*K/2*log(2*pi) + N*(K-J)/2*log(s2) + N/2*sum(log(delta(1:J))) + N*J/2 + N/2/s2*sum(delta(J+1:end));
detail = L + params/2;

%% codelength
msglen = assertion + detail;

end
