function [ result, combdes ] = CCVNSAlg( desupper, dlow, maxiter )

% CCVNSAlg: Variable Neighborhood Search Algorithm based on Changing Columns. 

% Description:

% The CC/VNS Algorithm is the local search algorithm and consists of two movements
% 1: Sign-switching any column
% 2: Swap any two columns and sign-switch one of them

% The four neighborhood used are:
% 1: Sign-switch any column
% 2: Swap any two columns
% 3: Sign-switch any two columns 
% 4: Choose three items and shift them to the right

% Objective: Find the best column structure for the lower design that
% minimizes a linear combination between the sum of squared correlation among
% two-factor interaction columns and between the pure-quadratic columns
% and the two-factor interactions

% Input Parameters:
% Cmatup = conference matrix used as upper matrix
% Cmatlow = conference matrix used as lower matrix
% a: parameter to denote the preference of the user (between 0 and 1)
% maxiter = number of iterations

% DATE: 21-SEP-2017

[n, m] = size(desupper); % Number of factors and runs
N = 2*n; % Run size of basic cOMARS design

% Allocate objects for the algorithm. Increase performance.
% Data for the neighborhoods

fact_vec = 1:m;
None = nchoosek(m, 2);
Ntwo = nchoosek(m, 3);
NSelTwo = nchoosek(fact_vec, 2);

% Variables to save the results
result = zeros(1, maxiter); % Objective value
combdes = cell(maxiter,1);

parfor jj = 1:maxiter
    
    % Extra for parallel computations
    maindeslower = dlow;
    NSelThree = nchoosek(fact_vec, 3);
    
    % Step 1: Starting Design
    krand = randsample(1:m, 1); 
    randchange = randsample(1:m, krand); % Random change of 'krand' columns
    maindeslower(:, randchange) = -1*maindeslower(:, randchange); % Sign-switch 'krand' columns
    selectperm = randsample(1:m, m); % Random permutation of the lower design
    startdeslower = maindeslower(:, selectperm); % Starting random lower design
    
    [objvalue, deslower] = CCAlgAfour(desupper, startdeslower, m, NSelTwo, None, N); % Apply CC algorithm and compute the first (current) objective value
    k = 1; % Initialize the Neighborhoods

    % Start VNS Algorithm
    while k <= 4 % Try to modify a little bit your actual solutions to see if there is a better solution in the nighbourhood

        % Select a neighborhood structure
        switch k
            
            %==============================================================
            case 1 % N1. Sign-switch any column 
                
                selectcolsNOne = randperm(m); % Random Order of the elements in the neighborhood
                k = 2; % If no improvement is found, move to the second neighborhood
                
                for ll = 1:m
                    
                    deslower(:, selectcolsNOne(ll)) = -1*deslower(:, selectcolsNOne(ll)); % Sign switch column 'll'

                    % Apply CC Algorithm
                    [resvalue, DesLowNeigh]  = CCAlgAfour(desupper, deslower, m, NSelTwo, None, N);

                    if resvalue < objvalue % Evaluate if an improvement has found
                       objvalue = resvalue; % Update Objective Value
                       deslower = DesLowNeigh; % Update Lower Desiggn
                       k = 1; % Go back to the first neighborhood
                       break % First improvement algorithm
                    end

                    % Re-store the current design
                    deslower(:, selectcolsNOne(ll)) = -1*deslower(:, selectcolsNOne(ll));      
                end
                
                
            %==============================================================
            case 2 % N2. Swap any two columns 
                
                selectcolsNtwo = randperm(None); % Random Order of the elements in the neighborhood
                k = 3; % If no improvement is found, move to the next neighborhood structure
                for ll = 1:None
                    
                    swaptwocols = NSelTwo(selectcolsNtwo(ll), :); % Select the first element in the neighborhood
                    neworder = fact_vec;
                    neworder(swaptwocols) =  neworder(flip(swaptwocols));
                    
                    % Apply CC algorithm
                    [resvalue, DesLowNeigh]  = CCAlgAfour(desupper, deslower(:, neworder), m, NSelTwo, None, N);

                    if resvalue < objvalue
                       objvalue = resvalue; % Update Objective Value
                       deslower = DesLowNeigh; % Update Lower Desiggn
                       k = 1; % Go back to the first neighborhood
                       break % First improvement algorithm
                    end
                          
                end
               
            %==============================================================
            case 3 % N3. Sign-switch any two columns 
                
                selectcolsNtwo = randperm(None); % Random Order of the elements in the neighborhood
                k = 4; % If no improvement is found, move to the next neighborhood structure

                for ll = 1:None
                    
                    swaptwocols = NSelTwo(selectcolsNtwo(ll), :); % Select the first element in the neighborhood
                    deslower(:, swaptwocols) = -1*deslower(:, swaptwocols); % Sign switch of two columns

                    % Apply CC algorithm
                    [resvalue, DesLowNeigh]  = CCAlgAfour(desupper, deslower, m, NSelTwo, None, N);

                    if resvalue < objvalue
                       objvalue = resvalue; % Update Objective Value
                       deslower = DesLowNeigh; % Update Lower Desiggn
                       k = 1; % Go back to the first neighborhood
                       break % First improvement algorithm
                    end
                    % Re-store design
                    deslower(:, swaptwocols) = -1*deslower(:, swaptwocols); % Sign switch of two columns
                          
                end
                
                
             %=============================================================
             case 4 % N4. Choose three items and shift them to the right
                
                selectcolsNthree = randperm(Ntwo); % Random order of the elements in the neighborhood
                k = 5; % If no improvement is found, move to the next neighborhood structure

                for ll = 1:Ntwo

                    selthreecols = NSelThree(selectcolsNthree(ll), :);
                    shifttotheright = selthreecols([3 1 2]); % Shift elements one position  to the right
                    neworder = fact_vec;
                    neworder(selthreecols) =  neworder(shifttotheright);
                    
                    % Apply CC algorithm
                    [resvalue, DesLowNeigh]  = CCAlgAfour(desupper, deslower(:, neworder), m, NSelTwo, None, N);

                    if resvalue < objvalue
                       objvalue = resvalue; % Update Objective Value
                       deslower = DesLowNeigh; % Update Lower Desiggn
                       k = 1; % Go back to the first neighborhood
                       break % First improvement algorithm
                    end  

                end
                
        end % end of Switch


    end

    % Save the best combined design
    % Add the center points
    %DD = [desupper; -1*desupper; deslower; -1*deslower; zeros(nzero, m)];
    %result(jj) = lcSOcorr( DD, NSelTwo, None, size(DD,1));
    result(jj) = lcSOcorr( [desupper; deslower], NSelTwo, None, N);
    combdes{jj} = [desupper; deslower];
end


end