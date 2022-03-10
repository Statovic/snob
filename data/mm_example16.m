% Example of mixture models of geometric distributions
clear;

%% Generate data
N = 1e2;    % number of data points
K = 3;      % number of classes
w = [0.2, 0.5, 0.3]; % mixing proportions
w = w ./ sum(w);     % make sure these sum up to 1
p = [0.7, 0.2, 0.4]; % parameter of each geometric distribution

x = nan(N, 1);
TrueInd = nan(N,1);
counts = mnrnd(N, w);
ix_start = 1;
for i = 1:K
    ix_end = ix_start + counts(i) - 1;
    x(ix_start:ix_end) = geornd(p(i), counts(i), 1);
    TrueInd(ix_start:ix_end) = i;
    ix_start = ix_end + 1;
end

%% Fit mixture model
mm = snob(x, {'geometric',1}, 'k', 3);
mm_PlotModel1d(mm, x, 1);

% Print Rand index and Adjusted Rand Index (ARI)
[~,EstInd] = max(mm.r, [], 2);          % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd); % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')