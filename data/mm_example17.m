%% Example - Simple mixture of regression models
clear;

%% Generate data
rng("default");
n = 50;
x1 = randn(n,1) * 1;
x2 = randn(n,1) * 5;
y1 = +5 * x1 + 1 + randn(n,1) * 3;
y2 = -3 * x2 + 5 + randn(n,1) * 5;
x = [x1; x2];
y = [y1; y2];

%% Plot data
subplot(1,3,1); plot(x1, y1,'rx'); xlabel('x1'); ylabel('y1'); grid;
subplot(1,3,2); plot(x2, y2, 'bx'); xlabel('x2'); ylabel('y2'); grid;
subplot(1,3,3); plot(x, y, 'kx'); xlabel('x'); ylabel('y'); grid;
hold on;

%% Fit mixture model
%  x ~ Normal(mu_1, sigma_1^2)
%  y ~ Normal(mu_2, sigma_2^2), mu2 = b0 + b1*x
mm = snob([x y], {'norm', 1, 'linreg',[2 1]}, 'k', 1);
mm_Summary(mm);

%% Plot mixture of regressions
N = length(y);
mu = zeros(N,mm.nClasses);
for k = 1:mm.nClasses
    m = mm.class{k}.model{end};    
    
    b0  = m.theta(2);
    b   = m.theta(3:end);
    mu(:,k)  = b0 + x*b;
end
plot(x,mu,'-','LineWidth',1.5);
hold off;

%% Compute adjusted rand index for the best model
TrueInd = [ones(n,1); ones(n,1)*2];
[~,EstInd] = max(mm.r, [], 2);              % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd);     % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);

mclass = minmis(mm, TrueInd);
fprintf('*** Minimum number of misclassifications is %d.\n', mclass);
fprintf('\n')
