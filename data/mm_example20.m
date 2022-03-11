%% Example - Analyis of the aphids data from:
% Turner, T. R. (2000). 
% Estimating the Propagation Rate of a Viral Infection of Potato Plants via Mixtures of Regressions. 
% Journal of the Royal Statistical Society. Series C (Applied Statistics), 49(3), 371–384. 
% http://www.jstor.org/stable/2680771

% Start with a clean workspace
clear;

% Load the data set
load aphids.mat;
NumAphids = aphids(:,1);
NumInfected = aphids(:,2);

% Visualise the data
plot(NumAphids, NumInfected, 'ko', ...
    'MarkerSize',4, 'MarkerFaceColor','k','MarkerEdgeColor','k');
xlim([0 340]); ylim([-4 29]); grid on;
xlabel('Number of aphdis released in the flight chamber');
ylabel('Number of infected plants');
set(gcf,'color','w');
hold on;

% Fit a mixture of linear regression models
% Following Turner (2000), we do not model the first column (Num. aphids)
% and instead assume the second column (Num. infected) follows a 
% mixture of Gaussian linear regressions. We use the keyword 'skip' to
% tell snob which columns can be ignored. Alternatively, we may choose to
% model the first column as a, say, normal distribution.
% Note that we start with 4 classes and let snob automatically choose the 
% number of classes.
mm = snob(aphids, {'skip', 1, 'linreg', [2 1]}, 'k', 4, 'display', false);

% Plot fitted regression lines
N = 100;
n = length(NumAphids);
x = linspace(1,350,N)';
mu = zeros(N, mm.nClasses);
res = zeros(n, mm.nClasses);
fit = zeros(n, mm.nClasses);
for k = 1:mm.nClasses
    m = mm.class{k}.model{end};
    b0  = m.theta(2);
    b   = m.theta(3:end);
    mu(:,k)  = b0 + x*b;
    fit(:,k) = b0 + NumAphids*b;
    res(:,k) = (NumInfected - fit(:,k));
end
plot(x, mu, '--', 'LineWidth', 1, 'Color','b');
hold off;

% Determine size of each plot point
scale = mm.r;

% Diagnostic plot: number of aphids vs residuals
figure;
scatter(NumAphids, res, scale*10, 'k', 'filled');
xlim([0 340]); ylim([-22 24]); grid;
xlabel('Number of aphdis released in the flight chamber');
ylabel('Residuals');
grid on; box on;
set(gcf,'color','w');

% Diagnostic plot: fitted values vs residuals
figure;
scatter(fit, res, scale*10, 'k', 'filled');
xlim([0 27]); ylim([-8 24]); grid;
xlabel('Fitted values');
ylabel('Residuals');
grid on; box on;
set(gcf,'color','w');



