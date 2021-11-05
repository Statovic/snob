function mm_KLstats(mm, data)

nClasses = mm.nClasses;
if(nClasses > 1)
    % Fit a mixture model with one class for the entire sample
    Pop = snob(data, mm.opts.ModelList, 'k', 1, 'fixedstructure', true, 'display', false);

    % For each attribute
    nModels = mm.nModelTypes;
    for i = 1:nModels
        klmat = zeros(nClasses+1, nClasses+1);
        
        % Compute KL matrix between class j and k for attr i 
        for j = 0:nClasses
            for k = 0:nClasses
                if(j ~= k)
                    klmat(j+1,k+1) = ComputeKL(Pop, mm, i, j, k);
                end
            end
        end

        % Print KL matrix
        fprintf('** %s ~ %s distribution\n', mm.opts.VarNames{i}, mm.ModelTypes{i}.type);
        printMat(klmat);
        fprintf('\n');
    end

end


end

%% Compute KL between two models
function kl = ComputeKL(Pop, mm, attr, classSrc, classTgt)

% Attribute type
type = mm.ModelTypes{attr}.type;

% Parameters of source and target model
thetaSrc = Pop.class{1}.model{attr}.theta;
if(classSrc > 0)
    thetaSrc = mm.class{classSrc}.model{attr}.theta;
end
thetaTgt = Pop.class{1}.model{attr}.theta;
if(classTgt > 0)
    thetaTgt = mm.class{classTgt}.model{attr}.theta;
end

% Compute KL
switch type
    case 'exp'
        p0 = thetaSrc; p1 = thetaTgt;
        kl = -1 + p0/p1 - log(p0) + log(p1);
    
    case 'Gaussian'
        m0 = thetaSrc(1); v0 = thetaSrc(2);
        m1 = thetaTgt(1); v1 = thetaTgt(2);
        kl = 0.5*(v0/v1 + (m1 - m0)^2/v1 + log(v1/v0) - 1);

    case 'multi'
        kl = sum(thetaSrc .* log(thetaSrc ./ thetaTgt));

    case 'mvg'
        d = mm.ModelTypes{attr}.nDim;
        m0 = thetaSrc(1:d);
        S0 = reshape(thetaSrc(d+1:end),d,d);
        m1 = thetaTgt(1:d);
        S1 = reshape(thetaTgt(d+1:end),d,d);        

        kl = mvgkl(m0, S0, m1, S1) / d; % per-dimension

    case 'Poisson'
        p0 = thetaSrc; p1 = thetaTgt;
        kl = -p0 + p1 + p0*log(p0) - p0*log(p1);

    otherwise
        error('KL divergence not implemented for type of model');
end


end

%% Print the KL matrix
function printMat(x)

n = size(x,1);

fprintf('%15s','');
for i = 1:n
    str = 'Pop.';
    if(i>1)
        str = ['Class ',num2str(i-1)];
    end
    fprintf('%15s', str);
end
fprintf('\n');

for i = 1:n
    str = 'Pop.';
    if(i>1)
        str = ['Class ',num2str(i-1)];
    end
    fprintf('%15s', str);
    for j = 1:n
        if(i~=j)
            fprintf('%15.4f', x(i,j));
        else
            fprintf('%15s', '-');
        end
    end
    fprintf('\n');
end


end

