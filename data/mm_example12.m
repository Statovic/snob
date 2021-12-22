%% Example - Simulated data (mixtures of multivariate Gaussian distributions)
clear;

rng(1);

%% Generate data
K = 3;          % number of classes
D = 2;          % dimension of each data point
n = ones(1,K)*25;   % sample size of each class

muMin = -5;    % Each mu \in [muMin, muMax]
muMax = +5;
muRange = muMax - muMin;

data = [];
Mu = cell(K,1);
Sigma = cell(K,1);
TrueInd = [];
for i = 1:K
    Sigma{i} = randcorr(D);
    Mu{i} = rand(D,1)*muRange + muMin;
    
    data  = [data; mvnrnd(Mu{i}, Sigma{i}, n(i))];
    TrueInd = [TrueInd; i*ones(n(i),1)];
end

%% Mixture model
mm = snob(data, {'mvg',1:D}, 'k', 1);
mm_Summary(mm);

%% Print matrix of KL divergences
mm_KLstats(mm, data);

% Print Rand index and Adjusted Rand Index (ARI)
[~,EstInd] = max(mm.r, [], 2);          % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd); % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')