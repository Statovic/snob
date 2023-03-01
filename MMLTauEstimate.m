% Computes MML estimates of tau (roots of a polynomial)
function s2_mml = MMLTauEstimate(N, K, J, delta)

W = delta(1:J);                 % Top J eigenvalues
R = sum( delta((J+1):end) );    % Remainder    

% Polynomial coefficients stored in x0
x0 = zeros(1,J+2);  

% We start by setting a0 and a_n
x0(1) = (N*(K-J) - J*(J-1)) * (-1)^J;
x0(end) = -N*prod(W)*R;

% We can compute the polynomial coefficients manually for J \in {1,2,3,4}
% For J > 4, we use out general formula in the paper
switch J
    % J = 1
    case 1 
        x0(2) = N*R + W*(K*(N - 1) - N);
    % J = 2
    case 2
        x0(2) = -N*R - sum(W)*(K*(N - 1) - 2*N - 1);
        x0(3) = N*R*sum(W) + W(1)*W(2)*(K*(N - 2) - 2*N);
    % J = 3
    case 3
        x0(2) = N*R + sum(W)*(K*(N - 1) - 3*N - 4);
        x0(3) = -N*R*sum(W) - (W(1)*W(2) + W(1)*W(3) + W(2)*W(3))*(K*(N - 2) - 3*N - 2);
        x0(4) = N*R*(W(1)*W(2) + W(1)*W(3) + W(2)*W(3)) + prod(W)*(K*(N - 3) - 3*N);
    % J = 4
    case 4
        x0(2) = -N*R - sum(W)*(K*(N - 1) - 4*N - 9);
        x0(3) = N*R*sum(W) + (W(1)*W(2) + W(1)*W(3) + W(1)*W(4) + W(2)*W(3) + W(2)*W(4) + W(3)*W(4))*(K*(N - 2) - 4*N - 6);
        x0(4) = -N*R*(W(1)*W(2) + W(1)*W(3) + W(1)*W(4) + W(2)*W(3) + W(2)*W(4) + W(3)*W(4)) - ...
            (W(1)*W(2)*W(3) + W(1)*W(2)*W(4) + W(1)*W(3)*W(4) + W(2)*W(3)*W(4))*(K*(N - 3) - 4*N  - 3);
        x0(5) = N*R*(W(1)*W(2)*W(3) + W(1)*W(2)*W(4) + W(1)*W(3)*W(4) + W(2)*W(3)*W(4)) + prod(W)*(K*(N - 4) - 4*N);
    otherwise
        ind = (J+1);
        
        %% General Formula [Below we have: j=1 linear term; j=2 quadratic term; j=3 cubic term, etc]
        for j = 1:J
            S1 = 1; S2 = 1;
            if((J-j) > 0)
                ix = nchoosek(1:J, J-j);
                S1 = sum( prod( reshape(W(ix),size(ix)), 2) );        
            end    
            if(J-j+1 > 0)
                ix = nchoosek(1:J, J-j+1);
                S2 = sum( prod( reshape(W(ix),size(ix)), 2) );
            end
        
            x0(ind) = (-1)^(j+1) * (N*R*S1 + S2*(K*(N-J+j-1) - J*N - (j-1)*(J-1)));
            ind = ind - 1;
        end        
end

%% Compute roots of polynomial with coefficients in x0
x0 = roots(x0);
s2_mml = MMLValidSolution(N, K, J, delta, x0);

end
