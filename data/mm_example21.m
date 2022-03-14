%% Example - Simulated data (mixtures of gamma distributions)
% Function 3 (Section 4.2) from 
% Wiper, M., Insua, D. R., & Ruggeri, F. (2001). 
% Mixtures of Gamma Distributions with Applications. 
% Journal of Computational and Graphical Statistics, 10(3), 440–454. 
% http://www.jstor.org/stable/1391098
clear;

% Seed the random number generator
rng("default");

% Generate n=400 data points from a mixture of three gamma distributions 
x = [gamrnd(40,1/20,80,1); gamrnd(6,1,240,1); gamrnd(200,1/20,80,1)];

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate gamma distribution: {'gamma',1}
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
%     (3) We start the search with 1 mixtures {'k',1}
mm = snob(x, {'gamma',1},'k',1);

% Print a summary of the discovered mixture model.
mm_Summary(mm);

% Plot the PDF of the mixture model
mm_PlotModel1d(mm, x, 1);

% Plot the CDF of the fitted mixture model 
minX = 0;    % range of values for plotting
maxX = 16;
nPts = 1e3;
x = linspace(minX, maxX, nPts)';
y = zeros(nPts, 1);

for k = 1:mm.nClasses
    prop = mm.a(k); % mixing proportion
    theta = mm.class{k}.model{1}.theta; % parameters for class k
    mu = theta(1);  % a = shape, b = scale
    phi = theta(2);
    a = phi; b = mu / phi;
    y = y + (prop * gamcdf(x,a,b));   % CDF
end

figure;
plot(x, y, 'k-');
grid;
xlabel('X', 'fontsize', 16);
ylabel('F(X)', 'fontsize', 16);