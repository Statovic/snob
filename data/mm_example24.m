%% Example - Simulated data (mixtures of PCAs)
clear;

rng(1);

%% Generate data
K = 3;          % number of classes
D = 8;          % dimension of each data point
J = 2;          % number of PCs 
n = [200 200 100 100];   % sample size of each class

muMin = -1;    % Each mu \in [muMin, muMax]
muMax = +1;
muRange = muMax - muMin;

data = [];
Mu = cell(K,1);
Sigma = cell(K,1);
Alpha = cell(K,1);
TrueInd = [];
for i = 1:K
    alpha = sort(abs(trnd(1,D,1)),'descend'); % Factor lengths
    s2 = mean(alpha((J+1):end));              % Residual variance sigma^2
    alpha = alpha(1:J);
    R = randn(D,J);     % Sample directions uniformly from a unit n-sphere
    R = bsxfun(@rdivide, R, sqrt(sum(R.^2)));
    A = R*diag(alpha);
    Alpha{i} = alpha;
    Mu{i} = rand(D,1)*muRange + muMin;
    Sigma{i} = A*A' + s2*eye(D);   % Variance-covariance matrix
    
    % Sample data
    data     = [data; mvnrnd(Mu{i}, Sigma{i}, n(i))];
    
    TrueInd = [TrueInd; i*ones(n(i),1)];
end

%% Mixture model
mm = snob(data, {'pca',1:D}, 'k', 1, 'useparallel', true, 'numpcs', J);
mm_Summary(mm);

%% Print matrix of KL divergences
mm_KLstats(mm, data);

%% Print Rand index and Adjusted Rand Index (ARI)
[~,EstInd] = max(mm.r, [], 2);          % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd); % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')
