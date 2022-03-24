%% Example - Simulated data (mixtures of Pareto Type II distributions)
clear;

rng('default');

%% Generate data
K = 2;               % number of classes
D = 2;
n = ones(1,K)*100;   % sample size of each class

data = [];
Mu = cell(K,1);
TrueInd = [];
for i = 1:K
    Mu{i} = abs(trnd(1,2,1));
    
    data  = [data; prtrnd(n(i),Mu{i}(1),1,Mu{i}(2))];
    TrueInd = [TrueInd; i*ones(n(i),1)];
end

%% Mixture model
mm = snob(data, {'pareto2',1}, 'k', 1);
mm_Summary(mm);

%% Print matrix of KL divergences
mm_KLstats(mm, data);

%% Print Rand index and Adjusted Rand Index (ARI)
[~,EstInd] = max(mm.r, [], 2);          % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd); % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')
