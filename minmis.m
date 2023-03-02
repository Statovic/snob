%MINMIS    Computes the minimum number of misclassifications
%  MINMIS(.) computes the minimum number of missclassifications by the
%  mixture model by label rotation given the true class labels. The total 
%  number of estimated classes (K) cannot exceed K=8 as there are K! 
%  labellings to check; e.g., 8! = 40,320.
%  
%  The input arguments are:
%   mm       - mixture model structure as returned by snob
%   TrueInd  - true class labels
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function minCl = minmis(mm, TrueInd)

% Estimate class assignments by max probability
[~,EstInd] = max(mm.r, [], 2);          
Kpred = length(unique(EstInd));

% Check how many classes exist
if(Kpred > 8)
    error(['Too many rotations to check [Max Rot = ', num2str(factorial(Kpred)), ']']);
end

%% Compute minimum number of misclassifications by label rotation
Labels = perms(1:Kpred);
m = length(Labels);

MisClass = zeros(m,1);
for i = 1:m
    lbl = Labels(i,:);
    EstLabel = lbl(EstInd);
    MisClass(i) = sum( EstLabel(:) ~= TrueInd(:) );
end

minCl = min(MisClass);

end
