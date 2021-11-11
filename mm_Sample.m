function [x, TrueInd] = mm_Sample(N, model, weights, varargin)

%% Input checking
if(any(mod(N,1)))
    error('#samples must be an integer greater than 0');
end
if(any(weights < 0))
    error('mixing proportions must be positive numbers')
end
params = cell2mat(varargin);
if(any(isinf(params(:))) || any(isnan(params(:))))
    error('Inf or NaN values detected in parameters')
end
if(length(weights) ~= length(varargin))
    error('incorrect number of parameters')
end

%% Check model
K = length(weights);    % number of classes
D = 1;                  % data dimension
for i = 1:K
    theta = varargin{i};
    k = length(varargin{1});
    m = length(theta);

    switch model
        case {'exp','geometric','poisson'}
            posCheck = true;               
    
        case {'beta','gamma','igauss','negb','weibull'}
            posCheck = [true, true];

        case 'multi'
            posCheck = true(k,1);

        case 'mvg'
            D = 1/2*(sqrt(1 + 4*k) - 1);            
            posCheck = false(k,1);

        case 'sfa'
            D = k / 3;
            posCheck = false(k,1);

        case {'laplace','norm'}
            posCheck = [false, true];
    
    otherwise
        error('model not implemented');        
    end    

    checkM = length(posCheck);
    if(m ~= checkM)
        error('model parameters not correctly specified: too many or too few');
    end
    if(any(theta(posCheck)<0))
        error('some model parameter must be positive')
    end
end


%% Generate data
x = nan(N,D);
TrueInd = nan(N,1);

weights = weights ./ sum(weights);  % make sure the weights sum up to 1
counts = mnrnd(N, weights);
ix_start = 1;
for i = 1:K
    ix_end = ix_start + counts(i) - 1;
    
    x(ix_start:ix_end, 1:D) = GenData(model, counts(i), varargin{i});
    TrueInd(ix_start:ix_end) = i;

    ix_start = ix_end + 1;
end

end

% Do the actual data generation
function x = GenData(model, N, params)

switch model
    case 'beta'
        a = params(1); b = params(2);
        x = betarnd(a, b, N, 1);
        x = min(max(eps,x), 1-eps);
    case 'exp'
        theta = params;
        x = exprnd(theta, N, 1);        
    case 'gamma'
        mu = params(1); phi = params(2);
        k = phi; theta = mu / phi;
        x = gamrnd(k, theta, N, 1);
    case 'geometric'
        theta = params;
        x = geornd(theta, N, 1);    
    case 'igauss'
        mu = params(1); lambda = params(2);
        x = randinvg(mu*ones(N,1), lambda);
    case 'laplace'
        mu = params(1); b = params(2);
        x = laprnd(mu, b, N);
    case 'multi'  % not required
        theta = params;
        counts = mnrnd(1, theta, N);
        x = ones(N,1);
        for j = 2:size(counts,2)
            ix = (counts(:,j) == 1);
            x(ix) = j;
        end
    case 'mvg'
        k = length(params);
        D = 1/2*(sqrt(1 + 4*k) - 1);   % D variate normal
        mu = params(1:D);
        Sigma = reshape(params(D+1:end),D,D);
        [~,posdef] = cholcov(Sigma,0);
        if(posdef)
            error('Sigma not positive definite');
        end
        x = mvnrnd(mu, Sigma, N);
    case 'negb'
        r = params(1); p = params(2);
        x = nbinrnd(r, p, N, 1);
    case 'norm'
        mu = params(1); v = params(2);
        x = normrnd(mu, sqrt(v), N, 1);
    case 'poisson'
        theta = params;
        x = poissrnd(theta, N, 1);       
    case 'sfa'        
        k = length(params);
        D = k/3;   % D variate "normal"
        mu = params(1:D);      
        sigma = params(D+1:2*D); 
        a_sfa = params(2*D+1:end); 
        Sigma = a_sfa*a_sfa' + diag(sigma.^2);    
        [~,posdef] = cholcov(Sigma,0);
        if(posdef)
            error('Sigma not positive definite');
        end        
        x = mvnrnd(mu, Sigma, N);
    case 'weibull'
        lambda = params(1); k = params(2);
        x = wblrnd(lambda, k, N, 1);       
end

end
