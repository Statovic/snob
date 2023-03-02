%% Example - Diabetes data
clear;

% Load the data. 
% The data consists of 10 predictors (X = SEX, AGE, BMI, BP and six blood serum measurements S1, S2, S3, S4, S5, S6) of diabetes (Y). 
load data/diabetes;
X = [SEX, AGE, BMI, BP, S1, S2, S3, S4, S5, S6];
[n,p] = size(X);
data = [X, Y];

% Run Snob with the following options: 
%     (1) the data is modelled as follows: {'multi',1,'mvg',2:4,'mvg',5:10,'linreg',[11,1:p]}
%               'multi',1           column 1 (SEX) follows a multinomial distribution
%               'mvg',2:4           columns 2,3,4 (AGE,BMI,BP) follow a multivariate Gaussian distribution
%               'mvg',5:10          columns 5:10 (S1-S6) follow a multivaraite Gaussian distribution
%               linreg',[11,1:p]    column 11 is a Gaussian linear regression; covariates are columns 1:10
%     (2) initial number of classes is 4 ('k',4)
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
VarNames = {'SEX', 'AGE', 'BMI', 'BP', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'DIABETES'};
mm = snob(data, {'multi',1,'mvg',2:4,'mvg',5:10,'linreg',[11,1:p]},'k',4,'varnames', VarNames);

% Print a summary of all the model parameters.
% There are two classes. The total message length is ~14,867 nits.
mm_Summary(mm);

% Print a matrix of KL divergences for the SFA model
% Label 'Pop.' denotes the single class model (ie, no mixture)
fprintf('KL divergences between each class:\n');
mm_KLstats(mm, data);

% Now, suppose we model AGE,BMI and BP as uncorrelated [spherical Gaussian]
mm2 = snob(data, {'multi',1,'norm',2:4,'mvg',5:10,'linreg',[11,1:p]}, ...
    'k',4,'varnames', VarNames, 'display', false);
mm_Summary(mm2);

% The codelengths are...
[mm.msglen, mm2.msglen]

% The uncorrelated model is 
exp( -(mm2.msglen - mm.msglen) )
% times more likely a posteriori than the MVG model.

% Show fitted values vs residuals
res = zeros(n, mm2.nClasses);
fit = zeros(n, mm2.nClasses);
nc = mm.Nk;
scale = mm2.r;
for k = 1:mm2.nClasses
    m = mm2.class{k}.model{end};
    b0  = m.theta(2);
    b   = m.theta(3:end);
    fit(:,k) = b0 + X*b;
    res(:,k) = (Y - fit(:,k));
end

figure;
scatter(fit, res, scale*5, 'k', 'filled');
grid;
xlabel('Fitted values');
ylabel('Residuals');
grid on; box on;
set(gcf,'color','w');