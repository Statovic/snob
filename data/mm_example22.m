%% Example - Simulated data (mixtures of dirichlet distributions)
clear;

rng('default');

%% Generate data
K = 2;          % number of classes
D = 4;          % dimension of each data point
n = ones(1,K)*150;   % sample size of each class

muMin = 1e-2;    % Each theta \in [muMin, muMax]
muMax = 1e+2;
muRange = muMax - muMin;

data = [];
Mu = cell(K,1);
TrueInd = [];
for i = 1:K
    Mu{i} = rand(D,1)*muRange + muMin;
    
    data  = [data; dirrnd(Mu{i}, n(i))];
    TrueInd = [TrueInd; i*ones(n(i),1)];
end

%% Mixture model
mm = snob(data, {'dirichlet',1:D}, 'k', 1);
mm_Summary(mm);

%% Print matrix of KL divergences
mm_KLstats(mm, data);

%% Print Rand index and Adjusted Rand Index (ARI)
[~,EstInd] = max(mm.r, [], 2);          % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd); % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')
