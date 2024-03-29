function mm = mm_EM(mm, data)

%% do search
done = false;
iter = 1;
L = zeros(mm.opts.emmaxiter, 1);

while not(done)
    
    %% Collapse mixtures if required ...
    if(min(mm.Nk) < mm.MinMembers)             
        [mm, removed] = mm_Collapse(mm, mm.Nk < mm.MinMembers);
        if(mm.opts.display && removed > 0)
            fprintf('             Removing %d class(es) due to insufficient membership...\n', removed);
        end           
    end
    
    %% Estimate class proportions
    [mm.a, mm.r, mm.Nk, ~] = mm_EstimateMixingWeights(mm, data);
        
    if(min(mm.Nk) >= mm.MinMembers || mm.nClasses == 1) % enough data to fit MM?
        %% Estimate models parameters
        mm = mm_EstimateTheta(mm, data, 1:mm.nClasses);

        %% New negative log-likelihood
        p = mm_Likelihood(mm, data, 1:mm.nModelTypes);
        p = max(min(p, 700), -700);     % log(p) \in [-700, 700]
        L(iter) = -sum(mylogsumexp(-p));  % OLD: -sum(log(sum(p,2)));
        
        %% Check exit conditions
        cond = (iter > 5) && mean( abs(diff(L(iter-4:iter))) ) < 1e-2;
        if iter >= mm.opts.emmaxiter || cond || any(isnan(mm.Nk))
            done = true;
        end       

        % next iteration
        iter = iter + 1;
    end
    
end

%% Get the final message length of the model
mm = mm_MsgLen(mm, data);

end
