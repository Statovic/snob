%% Example - Simulated data (mixtures of exponential distributions with fixed type I censoring)
clear;

% Seed the random number generator
rng(1);

%% We first generate some data with censoring.
% There are two classes:
C = 10;  % fixed censoring point
% (1) T ~ Exp(2), Y = min(T,C), Delta = I(T<=C)
n = 1e2;
T = exprnd(2, n, 1);
y1 = min(T, C);
delta1 = (T <= C)*1;    

% (2) T ~ Exp(1/5), Y = min(T,C), Delta = I(T<=C)
n = 50;
T = exprnd(1/5, n, 1);
y2 = min(T, C);
delta2 = (T <= C)*1;  

data = [[y1;y2], [delta1;delta2] ];

%% Run snob 
% We specify that we are dealing with censored exponential data.
% and start the search with k=1 classes.
mm = snob(data, {'cfixexp', [1,2]}, 'k', 1);

%% Lets look at the model snob discovered
mm_Summary(mm);

%% Example - Simulated data (mixtures of Weibull distributions with fixed type I censoring)
% There are two classes:
C = 2;  % fixed censoring point
n = 100;
T = wblrnd(1, 1, n, 1);
y1 = min(T, C);
delta1 = (T <= C)*1;    

n = 150;
T = wblrnd(2, 10, n, 1);
y2 = min(T, C);
delta2 = (T <= C)*1;  

data = [[y1;y2], [delta1;delta2] ];

%% Run snob 
% We specify that we are dealing with censored exponential data.
% and start the search with k=1 classes.
mm = snob(data, {'cfixweibull', [1,2]}, 'k', 5);

%% Lets look at the model snob discovered
mm_Summary(mm);
