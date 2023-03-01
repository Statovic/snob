%MINMIS    Computes the minimum number of misclassifications
%  MINMIS(.) computes the minimum number of missclassifications by the
%  mixture model given the true class labels.
%  
%  The input arguments are:
%   mm       - mixture model structure as returned by snob
%   TrueInd  - true class labels
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function minCl = minmis(mm, TrueInd)

% Estimated class assignments by max probability
[~,EstInd] = max(mm.r, [], 2);          

% compute min misclassification by label rotation
Kpred = length(unique(EstInd));
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
