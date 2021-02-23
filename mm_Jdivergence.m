function Jdiv = mm_Jdivergence(mm, wModels)

% If model matrix wasn't passed
if(~exist('wModels','var'))
    wModels = 1:mm.nModelTypes;
end

M = length(wModels);
K = mm.nClasses;   % number of mixtures
Jdiv = cell(M,1);

%% Compute J divergence 
% For each model...
for i = wModels
    
    Jdiv{i} = zeros(K,K);
    
    % Compute J divergence matrix between all the classes
    for j = 1:(K-1)
        for k = (j+1):K
            
            switch mm.ModelTypes{i}.type 
                % Exponential distribution
                case 'exp'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameter
                    b = mm.class{k}.model{i}.theta; % Model 2 parameter
                    
                    f = 1/2*(-2 + a/b + b/a);   % J divergence
                    
                % Geometric distribution
                case 'geometric'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameter
                    b = mm.class{k}.model{i}.theta; % Model 2 parameter
                    
                    f = 1/2*-(((a - b)*(log(1 - a) - log(1 - b)))/(a*b));
                    
                % Poisson distribution    
                case 'Poisson'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameter
                    b = mm.class{k}.model{i}.theta; % Model 2 parameter
                    
                    f = 1/2*(a - b)*log(a/b);
                    
                % Normal distribution
                case 'Gaussian'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    mu1 = a(1); v1 = a(2);
                    mu2 = b(1); v2 = b(2);
                    
                    f = (mu1-mu2)^2/2*(1/v1+1/v2) + 1/2*(v1/v2 + v2/v1 - 2);
                    
                % Weibull distribution
                case 'weibull'
                    a = mm.class{j}.model{i}.theta; % Model 1 parameters
                    b = mm.class{k}.model{i}.theta; % Model 2 parameters
                    lam1 = a(1); k1 = a(2);
                    lam2 = b(1); k2 = b(2);
                    
                    g = -psi(1);    % Euler–Mascheroni constant
                    f = log(k1) - k1*log(lam1) - log(k2) + k2*log(lam2) + (k1-k2)*(log(lam1) - g/k1) + (lam1/lam2)^k2*gamma(k2/k1+1) - 1;
                    f = f + log(k2) - k2*log(lam2) - log(k1) + k1*log(lam1) + (k2-k1)*(log(lam2) - g/k2) + (lam2/lam1)^k1*gamma(k1/k2+1) - 1;
                    f = f / 2;
                
                otherwise
                    error('Model not available');                    
            end
            
            % Save result to matrix
            Jdiv{i}(j,k) = f;
            Jdiv{i}(k,j) = f;
        end
    end
end


end