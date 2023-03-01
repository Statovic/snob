function [alpha, R, s2, collapsed, totalVar] = mmlpca(x, N, J)

% special case if J = 0
if(J == 0)
    alpha = 0;
    R = [];
    s2 = std(x(:),[],1);
    collapsed = true;
    return;
end

% Sx eigenvectors and eigenvalues
K = size(x,2);
Sx = (x'*x)/N;                  % sample variance-covariance matrix
[U,S,~] = svd(Sx,0);            % SVD of Sx
delta = diag(S);                % singular values of Sx
W = delta(1:J);
totalVar = sum(delta);
          
R = U(:,1:J);   % factor orientations

%% MML estimate of sigma2
s2 = MMLTauEstimate(N, K, J, delta);
alpha = max(0, sqrt(W - s2));    

collapsed = false;
if(min(alpha) < 1e-3)
    % MML collapsed a PC; try with J-1 PCs
    collapsed = true;

    estFound = false;
    for j = (J-1):-1:1
        W = delta(1:j);         
        R = U(:,1:j);   % factor orientations
        
        s2 = MMLTauEstimate(N, K, j, delta);
        alpha = max(0, sqrt(W - s2));    

        if(min(alpha) > 1e-3)
            estFound = true;
            break;
        end

    end
    if(~estFound) % J = 0 model
        alpha = 0;
        s2 = sum( x(:).^2 ) / (N*K);
        R = [];
    end
end

end
