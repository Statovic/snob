function [data, TrueInd] = mm_SampleMixture(N, model)

%% Input checking
if(any(mod(N,1)))
    error('#samples must be an integer greater than 0');
end

%% Generate data
data = nan(N, model.Ncols);
TrueInd = nan(N,1);

wClass = mnrnd(N, model.a);    % decide which class to sample from
ix_start = 1;
for i = 1:model.nClasses
    Nclass = wClass(i); % how many samples from this class
    if(Nclass > 0)
        ix_end = ix_start + Nclass - 1;    
    
        data(ix_start:ix_end, :) = SampleClass(Nclass, model.Ncols, model.class{i}.model);
        TrueInd(ix_start:ix_end) = i;
    
        ix_start = ix_end + 1;    
    end
end

% Permute the data
permind = randperm(N);
data = data(permind,:);
TrueInd = TrueInd(permind);

end

function data = SampleClass(N, K, models)

data = nan(N, K);
for i = 1:length(models)
    mdl = models{i};
    theta = mdl.theta;

    switch mdl.type
        case 'beta'
            a = theta(1); b = theta(2);
            x = betarnd(a, b, N, 1);
            x = min(max(eps,x), 1-eps);
            data(:,mdl.Ivar) = x;            

        case 'crndexp'
            a = theta(1); b = theta(2);
            T = exprnd(b, N, 1);    
            C = exprnd(a, N, 1);   
            data(:,mdl.Ivar) = [min(T,C), (T <= C)*1];

        case 'dirichlet'
            data(:,mdl.Ivar) = dirrnd(theta, N);  

        case 'exp'
            data(:,mdl.Ivar) = exprnd(theta, N, 1);                

        case 'gamma'
            mu = theta(1); phi = theta(2);
            shape = phi; scale = mu / phi;
            data(:,mdl.Ivar) = gamrnd(shape, scale, N, 1);    

        case 'Gaussian'
            mu = theta(1); v = theta(2);
            data(:,mdl.Ivar) = normrnd(mu, sqrt(v), N, 1);    

        case 'geometric'
            data(:,mdl.Ivar) = geornd(theta, N, 1);   

        case 'invGaussian'
            mu = theta(1)*ones(N,1); lambda = theta(2);
            data(:,mdl.Ivar) = randinvg(mu, lambda);  

        case 'laplace'
            mu = theta(1); b = theta(2);
            data(:,mdl.Ivar) = laprnd(mu, b, N); 

        case 'lognorm'
            mu = theta(1); s2 = theta(2);
            data(:,mdl.Ivar) = exp( normrnd(mu, sqrt(s2), N, 1) );

        case 'multi'
            counts = mnrnd(1, theta, N);
            x = ones(N,1);
            for j = 2:size(counts,2)
                ix = (counts(:,j) == 1);
                x(ix) = j;
            end
            data(:,mdl.Ivar) = x;       

        case 'mvg'
            k = length(theta);
            D = 1/2*(sqrt(1 + 4*k) - 1);   % D variate normal
            mu = theta(1:D);
            Sigma = reshape(theta(D+1:end),D,D);
            [~,posdef] = cholcov(Sigma,0);
            if(posdef)
                error('Sigma not positive definite');
            end
            data(:,mdl.Ivar) = mvnrnd(mu, Sigma, N);

        case 'negb'
            mu = theta(1); phi = theta(2);
            r = phi; 
            p = 1-mu/(mu+phi);
            data(:,mdl.Ivar) = nbinrnd(r, p, N, 1);     

        case 'pareto2'
            sigma = theta(1); alpha = theta(2);
            data(:,mdl.Ivar) = prtrnd(N, sigma, 1, alpha);

        case 'poisson'
            lambda = theta(1);
            data(:,mdl.Ivar) = poissrnd(lambda, N, 1);     

        case 'sfa'
           d = length(mdl.Ivar);
            mu = theta(1:d);
            sigma = theta((d+1):(2*d));
            a_sfa = theta((2*d+1):end);
            Sigma = a_sfa*a_sfa' + diag(sigma.^2);
            data(:,mdl.Ivar) = mvnrnd(mu, Sigma, N);

        case 'vmf'
            kappa = theta(1); mu = theta(2:end);
            data(:,mdl.Ivar) = vmfrnd(kappa,mu,N);

        case 'weibull'
            lambda = theta(1); k = theta(2);
            data(:,mdl.Ivar) = wblrnd(lambda, k, N, 1);      

        otherwise
            error('Model not supported')
    end
end

end