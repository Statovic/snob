%% Example - Simulated data (mixtures of linear regressions)
clear;

%% Generate data from K mixtures of linear regression models
% Each linear regression model has D covariates
K = 4;          % number of classes
D = 2;          % number of covariates
DataSizes = [100, 200, 250, 500]; % possible data sizes for each class
ind = mnrnd(length(DataSizes)-1, ones(1,K)./K) + 1;
n = DataSizes(ind);   % actual data sizes 

bMin = -10;    % Each mu \in [muMin, muMax]
bMax = +10;
bRange = bMax - bMin;
snr = 5;

data = [];
TrueInd = [];
fprintf('Data generated from the mixture of %d regression models:\n', K);
for i = 1:K
    % Generate X    
    S = 1;
    if(D > 1)
        S = randcorr(D);
    end        
    X = mvnrnd(zeros(1,D), S, n(i));
    
    % Generate betas
    b0 = rand(1)*bRange + bMin;
    b = rand(D,1)*bRange + bMin;
    
    % Generate targets
    mu = b0 + X*b;
    s = max(0.1, var(mu) / snr);
    y = mu + randn(n(i),1)*s;

    fprintf('  y = %+5.2f', b0);
    for j = 1:D
        if(b(j) > 0)
            fprintf(' + ');
        else
            fprintf(' - ');
        end
        fprintf('%4.2f X%d', abs(b(j)), j);
    end
    fprintf(' + N(0, %3.2f^2)', s);
    fprintf('\t(n = %3d)', n(i));
    fprintf('\n');

    data  = [data; [X y]];
    TrueInd = [TrueInd; i*ones(n(i),1)];
end
fprintf('\n');


%% Fit mixture model
if(D == 1)
    mm = snob(data, {'norm', 1, 'linreg', [2 1]}, 'k', 1, 'display', false);    % 1 covariate only
else
    mm = snob(data, {'mvg', 1:D, 'linreg', (D+1):-1:1}, 'k', 1, 'display', false);    
end
mm_Summary(mm);

%% Compute adjusted rand index
[~,EstInd] = max(mm.r, [], 2);              % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd);     % Compute Rand index and ARI

%% Print KL stats
if(mm.nClasses > 1)
    fprintf('KL divergence matrix:\n')
    mm_KLstats(mm, data);
end

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')