%MM_SUMMARY    Print a summary of the model structure and parameters.
%  MM_SUMMARY(.) prints a summary of the mixture model structure (number of
%  classes; mixing proportions) and all the corresponding model parameters.
%  
%  The input arguments are:
%   mm    - structure respresenting the complete mixture model
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function mm_Summary(mm)

%% How many classes and models per class
nClasses = mm.nClasses;
nModels  = mm.nModelTypes;

fprintf('Minimum Message Length Mixture Model\n');
fprintf('Message Length  = %8.2f\n', mm.msglen);    
fprintf('MML Assertion   = %8.2f\n', mm.Atheta + mm.constant + mm.Ak + mm.Aa);      
fprintf('Neg. LogLike    = %8.2f\n', mm.L);    
fprintf('BIC             = %8.2f\n', mm.BIC);    
fprintf('AIC             = %8.2f\n', mm.AIC);    
fprintf('Classes         = %8d\n', mm.nClasses);    
fprintf('Parameters      = %8d\n', round(mm.nParams));    
fprintf('Observations    = %8d\n', mm.N);    
fprintf('\n');

%% Print class information
% For each class
for k = 1:nClasses
    
    fprintf('%5s%3d: [%10s = %4.1f%%;%3s = %6.1f]\n', 'Class', k, 'Proportion', mm.a(k)*100, 'Nk', mm.Nk(k));
    % For each model (across columns), get likelihoods
    for i = 1:nModels
        fprintf('%10s%3d', 'Model', i);
        
        model = mm.class{k}.model{i};           % model                
        theta = mm.class{k}.model{i}.theta;
        Ivar  = mm.class{k}.model{i}.Ivar;
        switch model.type
            
            %% skip
            case 'skip'
                fprintf('%20s', 'No model');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});                    

            %% Dirichlet distribution
            case 'dirichlet'
                fprintf('%20s', 'Dirichlet');
                fprintf('%10s ', mm.opts.VarNames{Ivar(1)});                
                d = length(theta);
                for j = 2:d
                    fprintf('%s ', mm.opts.VarNames{Ivar(j)});                
                end
                fprintf('\n');                
                fprintf('%43s = [', 'theta');     
                for j = 1:d
                    fprintf('%4.2f', theta(j));
                    if(j ~= length(theta))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                
                
            %% von Mises-Fisher distribution
            case 'vmf'
                fprintf('%20s', 'von Mises-Fisher');
                fprintf('%10s ', mm.opts.VarNames{Ivar(1)});                
                d = length(theta) - 1;
                for j = 2:d
                    fprintf('%s ', mm.opts.VarNames{Ivar(j)});                
                end
                fprintf('\n');
             
                fprintf('%43s = %7.2f\n', 'kappa', theta(1));
                fprintf('%43s = [', 'mu');                
                mu = theta(2:end);               
                % mu
                for j = 1:length(mu)
                    fprintf('%4.2f', mu(j));
                    if(j ~= length(mu))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                
                                

            %% beta distribution
            case 'beta'
                fprintf('%20s', 'beta');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f', 'a',theta(1),'b',theta(2));            
            
            %% Weibull distribution
            case 'weibull'
                fprintf('%20s', 'Weibull');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f', 'lambda',theta(1),'k',theta(2));
                
            %% Weibull distribution with fixed type I censoring
            case 'cfixweibull'
                fprintf('%20s', 'Weibull Type I');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f%10s = %7.2f', 'lambda',theta(1),'k',theta(2),'Cens.',mm.ModelTypes{i}.c);                
            
            %% Exponential with random type I censoring
            case 'crndexp'
                fprintf('%20s', 'Exp Rand');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f', 'alpha',theta(1),'beta',theta(2));
                
            %% Exponential with fixed type I censoring
            case 'cfixexp'
                fprintf('%20s', 'Exp Type I');
                fprintf('%10s %s\n', mm.opts.VarNames{Ivar(1)},mm.opts.VarNames{Ivar(2)});              
                fprintf('%43s = %7.2f%10s = %7.2f', 'theta',theta(1),'Cens.',mm.ModelTypes{i}.c);
                
            %% Exponential
            case 'exp'
                fprintf('%20s', 'Exp');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});              
                fprintf('%43s = %7.2f', 'lambda',theta(1));
            
            %% Multinomial
            case 'multi'                
                fprintf('%20s', 'Multinomial');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});               
                fprintf('%43s = [', 'p');
                for j = 1:length(theta)
                    fprintf('%4.2f', theta(j));
                    if(j ~= length(theta))
                        fprintf(' ');
                    end
                end
                fprintf(']');
            
            %% Gaussian
            case 'Gaussian'
                fprintf('%20s', 'Gaussian');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f', 'mu',theta(1),'std',sqrt(theta(2)));

            %% Gaussian
            case 'lognorm'
                fprintf('%20s', 'Lognormal');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f', 'mu',theta(1),'std',sqrt(theta(2)));                

            %% Laplace
            case 'Laplace'
                fprintf('%20s', 'Laplace');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});
                fprintf('%43s = %7.2f%10s = %7.2f', 'mu',theta(1),'b',theta(2));
                
            %% Gamma    
            case 'gamma'
                fprintf('%20s', 'Gamma');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});  
                fprintf('%43s = %7.2f%10s = %7.2f', 'mu',theta(1),'phi',theta(2));                
                
            %% PCA
            case 'pca'
                d = mm.ModelTypes{i}.nDim;
                numPCs = model.J;
                fprintf('%20s', ['PCA(k=',num2str(numPCs),')']);        
                fprintf('%10s ', mm.opts.VarNames{Ivar(1)});                
                for j = 2:d
                    fprintf('%s ', mm.opts.VarNames{Ivar(j)});                
                end
                fprintf('\n');     

                mu = theta(1:d);
                alpha_pca = theta((d+1):(d+numPCs));
                R_pca = reshape(theta((d+numPCs+1):(d+numPCs+d*numPCs)),d,numPCs);
                A_pca = R_pca * diag(alpha_pca);
                s2 = theta(end);

                Sigma = A_pca*A_pca' + diag(s2)*eye(d);
                s = sqrt(diag(Sigma));
                R = Sigma ./ (s*s');                

                % mu
                fprintf('%43s = [', 'mu');
                for j = 1:length(mu)
                    fprintf('%4.2f', mu(j));
                    if(j ~= length(mu))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                   

                % sigma
                fprintf('%43s = [', 'std');                
                for j = 1:length(s2)
                    fprintf('%4.2f', sqrt(s2(j)));
                    if(j ~= length(s2))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                       
                
                % a
                fprintf('%43s = [', 'a');                
                for j = 1:length(alpha_pca)
                    fprintf('%+5.2f', alpha_pca(j));
                    if(j ~= length(alpha_pca))
                        fprintf(' ');
                    end
                end    
                fprintf(']');         

                if(model.collapsed)
                    fprintf('*');
                else
                    fprintf(', %s = %.2f', '||a / std||', norm(alpha_pca./sqrt(s2)));  
                end
                fprintf('\n');

                % total variance explained
                fprintf('%43s = [', 'Var. Explained');                
                for j = 1:length(alpha_pca)
                    fprintf('%+5.2f', alpha_pca(j)^2/model.totalVar*100);
                    if(j ~= length(alpha_pca))
                        fprintf(' ');
                    end
                end    
                fprintf(']\n');                 

                % R
                fprintf('%43s = [', 'R');
                if(d <= 10)
                    for j = 1:d
                        if(j>1)
                            fprintf('%47s','[');
                        end
                        for jj = 1:d
                            fprintf('%+8.2f', R(j,jj));
                            if(jj < d)
                                fprintf(' ');
                            else
                                fprintf(']');
                            end
                        end
                        if(j<d)
                            fprintf('\n');
                        end
                    end                                
                else
                    fprintf(' ... ]\n');
                end
                
            %% Single factor analysis    
            case 'sfa'
                d = mm.ModelTypes{i}.nDim;
                fprintf('%20s', 'SingleFA');        
                fprintf('%10s ', mm.opts.VarNames{Ivar(1)});                
                for j = 2:d
                    fprintf('%s ', mm.opts.VarNames{Ivar(j)});                
                end
                fprintf('\n');                
                fprintf('%43s = [', 'mu');
                
                mu = theta(1:d);      
                sigma = theta(d+1:2*d); 
                a_sfa = theta(2*d+1:end); 
                Sigma = a_sfa*a_sfa' + diag(sigma.^2);
                s = sqrt(diag(Sigma));
                R = Sigma ./ (s*s');
                
                % mu
                for j = 1:length(mu)
                    fprintf('%4.2f', mu(j));
                    if(j ~= length(mu))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');       
                
                % sigma
                fprintf('%43s = [', 'std');                
                for j = 1:length(sigma)
                    fprintf('%4.2f', sigma(j));
                    if(j ~= length(sigma))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                       
                
                % a
                fprintf('%43s = [', 'a');                
                for j = 1:length(a_sfa)
                    fprintf('%+5.2f', a_sfa(j));
                    if(j ~= length(a_sfa))
                        fprintf(' ');
                    end
                end
                fprintf(']');                                 
                if(model.collapsed)
                    fprintf('*');
                else
                    fprintf(', %s = %.2f', '||a / std||', norm(a_sfa./sigma));  
                end
                fprintf('\n');
                
                % R
                fprintf('%43s = [', 'R');
                for j = 1:d
                    if(j>1)
                        fprintf('%47s','[');
                    end
                    for jj = 1:d
                        fprintf('%+8.2f', R(j,jj));
                        if(jj < d)
                            fprintf(' ');
                        else
                            fprintf(']');
                        end
                    end
                    if(j<d)
                        fprintf('\n');
                    end
                end
                
            
            %% Multivariate Gaussian
            case 'mvg'
                
                d = mm.ModelTypes{i}.nDim;
                fprintf('%20s', 'MVGaussian');
                fprintf('%10s ', mm.opts.VarNames{Ivar(1)});                
                for j = 2:d
                    fprintf('%s ', mm.opts.VarNames{Ivar(j)});                
                end
                fprintf('\n');
                fprintf('%43s = [', 'mu');
                
                mu = theta(1:d);
                Sigma = reshape(theta(d+1:end),d,d);                
                sigma = sqrt(diag(Sigma));
                R = Sigma ./ (sigma*sigma');
                
                % mu
                for j = 1:length(mu)
                    fprintf('%4.2f', mu(j));
                    if(j ~= length(mu))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                
                
                % sigma
                fprintf('%43s = [', 'std');
                for j = 1:length(sigma)
                    fprintf('%4.2f', sigma(j));
                    if(j ~= length(sigma))
                        fprintf(' ');
                    end
                end
                fprintf(']\n');                                
                
                % R
                fprintf('%43s = [', 'R');
                for j = 1:d
                    if(j>1)
                        fprintf('%47s','[');
                    end
                    for jj = 1:d
                        fprintf('%+8.2f', R(j,jj));
                        if(jj < d)
                            fprintf(' ');
                        else
                            fprintf(']');
                        end
                    end
                    if(j<d)
                        fprintf('\n');
                    end
                end
            
            %% Poisson    
            case 'Poisson'
                fprintf('%20s', 'Poisson');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});  
                fprintf('%43s = %7.2f', 'lambda',theta(1));
                
            %% Geometric    
            case 'geometric'
                fprintf('%20s', 'geometric');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});    
                fprintf('%43s = %7.2f', 'p',theta(1));                
            
            %% Inverse Gaussian    
            case 'pareto2'
                fprintf('%20s', 'Pareto (II)');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});  
                fprintf('%43s = %7.2f%10s = %7.2f', 'sigma',theta(1),'alpha',theta(2));            
                
            %% Inverse Gaussian    
            case 'invGaussian'
                fprintf('%20s', 'Inv-Gaussian');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});  
                fprintf('%43s = %7.2f%10s = %7.2f', 'mu',theta(1),'lambda',theta(2));
               
            %% Negative binomial    
            case 'negb'
                fprintf('%20s', 'Neg-binomial');
                fprintf('%10s\n', mm.opts.VarNames{Ivar(1)});  
                mu = theta(1); phi = theta(2);
                theta = [phi, 1-mu/(mu+phi)];
                fprintf('%43s = %7.2f%10s = %7.2f', 'r',theta(1),'p',theta(2));                
            
            %% Linear regression
            case 'linreg'
                p = mm.ModelTypes{i}.CovIx;
                fprintf('%20s', 'Linreg');
                fprintf('%10s ~ 1', mm.opts.VarNames{Ivar});              
                if(~isempty(p))
                    fprintf(' + ');
                end
                for j = 1:length(p)
                    fprintf('%s', mm.opts.VarNames{p(j)});                
                    if(j<length(p))
                        fprintf(' + ');
                    end
                end
                fprintf('\n');
                fprintf('%43s = %7.2f%10s = %7.2f\n', 'b0',theta(2),'std',sqrt(theta(1)));
                fprintf('%43s = [', 'b');

                for j = 1:length(p)
                    fprintf('%5.2f', theta(j+2));
                    if(j ~= length(p))
                        fprintf(' ');
                    end                    
                end
                fprintf(']');
            
            %% Logistic regression    
            case 'logreg'
                p = mm.ModelTypes{i}.CovIx;
                fprintf('%20s', 'Logreg');
                fprintf('%20s ~ 1', ['logOdds(', mm.opts.VarNames{Ivar}, ')']);              
                if(~isempty(p))
                    fprintf(' + ');
                end                
                for j = 1:length(p)
                    fprintf('%s', mm.opts.VarNames{p(j)});                
                    if(j<length(p))
                        fprintf(' + ');
                    end
                end        
                fprintf('\n');
                fprintf('%43s = %7.2f\n', 'b0',theta(1));
                fprintf('%43s = [', 'b');
                for j = 1:length(p)
                    fprintf('%5.2f', theta(j+1));
                    if(j ~= length(p))
                        fprintf(' ');
                    end                    
                end
                fprintf(']');                
                
            otherwise
                error('Unknown model type');
        end
        
        fprintf('\n');
    end
    fprintf('\n');
end

end