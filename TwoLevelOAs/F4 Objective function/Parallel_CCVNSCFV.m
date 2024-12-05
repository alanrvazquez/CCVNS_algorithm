function [ result, combdes] = Parallel_CCVNSCFV( desupper, pardeslower, maxiter, vecbigM)

% PARALLEL version of the CC/VNS algorithm for optimizing the F4 vector of 
% concatenated designs. 
%
% INPUTS:
% desupper      The N-by-m upper design with coded levels -1 and
%               +1 (matrix).
% maindeslower  The N-by-m lower design with coded levels -1 and 
%               +1 (matrix).
% maxiter       The number of iterations.
% vecbigM       The 1-by-nelF4 vector with the weight for each component of 
%               the F4 vector, where nelF4 is the number of elements in the 
%               F4 vector of a strength-3 design.
%
% OUTPUTS:
% result    A maxiter-by-1 vector with the f(D) values of the optimized 
%           concatenated design for each iteration, see Supplementary 
%           Section A.2. 
% combdes   A maxiter-by-1 cell array containing the 2N-by-(m+1) optimized
%           concatenated designs with an extra factor and coded levels
%           -1 and +1, for each iteration.
%
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

%====SET UP INPUTS AND INITIALIZE VARIABLES================================
[N, m] = size(desupper); 
twoproN = 2*N; 

% Allocate data for the neighborhoods.-------------------------------------
fact_vec = 1:m;
None = nchoosek(m, 2); 
Ntwo = nchoosek(m, 3);

% Allocate data for 2FI matrix computations.-------------------------------
Selmintwo = nchoosek(1:m, 2); 
mchoosetwo = nchoosek(m, 2);

% Variables for computing the frequencies of the J4-characteristics.-------
nelF4=ceil(twoproN/16); % The number of elements in the F4 vector of a 
                        % strength-3 design.
clind = nelF4:-1:1;
PossibJ4s = twoproN-(clind-1)*16; % Possible values of the J4-characteristics.

% Contribution of the upper design.----------------------------------------
Ndiag = diag(twoproN*ones(1, None));
UpInt = TwoFIMat(desupper, Selmintwo, mchoosetwo, N);
UTU = UpInt'*UpInt;

% Restrict N4 size for larger cases.---------------------------------------
if m > 24 
    Ntwo = 2024;
end

% Variables to save the results.-------------------------------------------
result = zeros(maxiter, 1);
combdes = cell(maxiter, 1);

parfor jj = 1:maxiter
    
    % Step 0: Allocate variables for parallel computations-----------------
    maindeslower = pardeslower;
    NSelTwo = nchoosek(fact_vec, 2);
    NSelThree = nchoosek(fact_vec, 3);
    
    % Step 1: Generate random plan for maindeslower.-----------------------
    krand = randsample(1:m, 1); 
    randchange = randsample(1:m, krand); 
    maindeslower(:, randchange) = -1*maindeslower(:, randchange); 
    selectperm = randsample(1:m, m); 
    
    startdeslower = maindeslower(:, selectperm); 
    
    % Step 2: Generate starting plan for maindeslower.---------------------
    %         Apply CC algorithm.
    [objvalue, deslower]  = CCAlgCFV(startdeslower, m, UTU, UpInt, Selmintwo, nelF4, vecbigM, Ndiag, PossibJ4s, None, N);
    
    % Step 3: Explore the neighborhood structures.-------------------------
    k = 1; 
    while k <= 4 

        % Select a neighborhood structure.
        switch k
            
            %============================================================== 
            case 1 % N1. Switch signs of any column.  

                selectcolsNOne = randperm(m); % Random Order of the elements 
                                              % in the neighborhood.
                k = 2; % If no improvement is found, move to N2.
                
                % Explore N1.
                for ll = 1:m
                    
                    % Construct neighboring solution.
                    lowdes = deslower; 
                    lowdes(:, selectcolsNOne(ll)) = -1*lowdes(:, selectcolsNOne(ll)); 

                    % Apply CC Algorithm.
                    [resvalue, DesLowNeigh] = CCAlgCFV(lowdes, m, UTU, UpInt, Selmintwo, nelF4, vecbigM, Ndiag, PossibJ4s, None, N);
                    
                    if resvalue < objvalue % If improvement. 
                       objvalue = resvalue; % Update best objective value.
                       deslower = DesLowNeigh; % Update best lower design.
                       k = 1; % Go back to N1.
                       break % Stop exploring N1. First improvement algorithm.
                    end
                          
                end
                
                
            %============================================================== 
            case 2 % N2. Swap any two columns 

                selectcolsNtwo = randperm(None); % Random Order of the 
                                                 % elements in the neighborhood.
                k = 3; % If no improvement is found, move to N3.
                
                % Explore N2.
                for ll = 1:None
                    % Construct neighboring solution. 
                    swaptwocols = NSelTwo(selectcolsNtwo(ll), :); 
                    neworder = fact_vec;
                    neworder(swaptwocols) =  neworder(flip(swaptwocols));
                    
                    % Apply CC algorithm.
                    [resvalue, DesLowNeigh] = CCAlgCFV(deslower(:, neworder), m, UTU, UpInt, Selmintwo, nelF4, vecbigM, Ndiag, PossibJ4s, None, N);
                    if resvalue < objvalue % If improvement.
                       objvalue = resvalue; % Update best objective value.
                       deslower = DesLowNeigh; % Update best lower design.
                       k = 1; % Go back to N1.
                       break % Stop exploring N2. First improvement algorithm.
                    end
                          
                end
               
            %============================================================== 
            case 3 % N3. Switch signs of any two columns.  

                selectcolsNtwo = randperm(None); % Random Order of the 
                                                 % elements in the neighborhood.
                k = 4; % If no improvement is found, move to N4.
                
                % Explore N3.
                for ll = 1:None
                    
                    % Construct neighboring solution.
                    swaptwocols = NSelTwo(selectcolsNtwo(ll), :); 
                    lowdes = deslower;
                    lowdes(:, swaptwocols) = -1*lowdes(:, swaptwocols); 

                    % Apply CC algorithm.
                    [resvalue, DesLowNeigh] = CCAlgCFV(lowdes, m, UTU, UpInt, Selmintwo, nelF4, vecbigM, Ndiag, PossibJ4s, None, N);
                    if resvalue < objvalue % If improvement.
                       objvalue = resvalue; % Update best objective value.
                       deslower = DesLowNeigh; % Update best lower design.
                       k = 1; % Go back to N1.
                       break % Stop exploring N3. First improvement algorithm.
                    end
                          
                end
                
                
             %=============================================================    
             case 4 % N4. Choose any subset of three columns, move the first
                    %     two columns one position to the right and move
                    %     the third column to position 1.

                selectcolsNthree = randperm(Ntwo); % Random order of the 
                                                   % elements in the neighborhood.
                k = 5; % If no improvement is found, exit the neighborhoods. 
                
                % Explore N4.
                for ll = 1:Ntwo
                    
                    % Construct neighboring solution.
                    selthreecols = NSelThree(selectcolsNthree(ll), :);
                    shifttotheright = selthreecols([3 1 2]); 
                    neworder = fact_vec;
                    neworder(selthreecols) =  neworder(shifttotheright);
                    
                    % Apply CC algorithm.
                    [resvalue, DesLowNeigh] = CCAlgCFV(deslower(:, neworder), m, UTU, UpInt, Selmintwo, nelF4, vecbigM, Ndiag, PossibJ4s, None, N);
                    
                    if resvalue < objvalue % If improvement.
                       objvalue = resvalue; % Update best objective value.
                       deslower = DesLowNeigh; % Update best lower design.
                       k = 1; % Go back to N1.
                       break % Stop exploring N4. First improvement algorithm.
                    end  

                end
                
        end % end of Switch.

    end % end of while.

    % Step 4. Save the best designs----------------------------------------
    finaldes = horzcat([-1*ones(twoproN/2, 1); ones(twoproN/2, 1)], [desupper; deslower]); % Add the extra factor. 
    result(jj, 1) = bigMFfour(finaldes, twoproN, vecbigM, Selmintwo, mchoosetwo); % Save the resulting best f(D) value.
    combdes{jj} = finaldes; % Save optimized design.
end

end