clear;

%% Load the old failthful data set
% The data consists of two numerical continuous variables.
%
% Reference: 
% Azzalini, A. and Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. Applied Statistics 39, 357–365.
load data/oldfaithful;

%% Cluster waiting times using snob
waiting = oldfaithful(:,2);
mm = snob(waiting, {'norm', 1},'k',2,'varnames',{'waiting'});

%% Summary
% Print a summary of all the components (parameters and structure) of the
% mixture model we have discovered. Snob discovered two classes
%
mm_Summary(mm);
mm_Table(mm);

% Plot the mixture distribution
mm_PlotModel1d(mm, waiting, 1);