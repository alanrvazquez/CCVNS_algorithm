function [objvalue, lowdes] = CCAlg(desupper, lowdes, m, ncombtwo, nm, N, vecbigM, nelLambda, Lambda)

% Column Change algorithm for lcSOcorr. 

% Description:
% The CCAlgorithm consists of two movements
% 1: Sign-switching any column
% 2: Swap any two columns and sign-switch one of them

% Objective: Find the best element for each position of the
% lower design

% Objective Function: lcSOcorr

% Input Parameters:
% desupper = Upper Design
% deslower = Lower Design
% m = number of factors

% DATE: 21-SEP-2017

% Compute current objective function

objvalue = lcSOcorr([desupper; lowdes], ncombtwo, nm, N, vecbigM, nelLambda, Lambda);

% Flags to stop the algorithm
thetazero=0;
bflag=0;

% Start local search algorithm
while bflag < 1
% Step 2 Improve
    for ii = 1:m

        % Sign-switch column 'ii'
        lowdes(:, ii) = -1*lowdes(:, ii); 
        % Evaluate design
        resvalue = lcSOcorr([desupper; lowdes], ncombtwo, nm, N, vecbigM, nelLambda, Lambda);

        if resvalue < objvalue % If an improvement is found, 
            % Update objective value
            objvalue = resvalue; % Update objective value
            % lowdes already updated 
            continue % Do not execute the following and go to 'ii+1'
        end
        % If fold-over does not work, restore lowdes
        lowdes(:, ii) = -1*lowdes(:, ii);
        
        % Try to fix the worst column starting from the best columns
        % Swap factor 'ii' and the rest
        st1 = ii+1;
        for jj = st1:m
   
            lowdes(:, [ii jj]) = lowdes(:, [jj ii]); % Swap two columns
            resvalueTwo = lcSOcorr([desupper; lowdes], ncombtwo, nm, N, vecbigM, nelLambda, Lambda); % Positive Version     
            lowdes(:, ii) = -1*lowdes(:, ii);
            resvalue = lcSOcorr([desupper; lowdes], ncombtwo, nm, N, vecbigM, nelLambda, Lambda); % Negative Version

            % Determine which one is the best change
            RandSwitch = 1;
            if resvalueTwo < resvalue                  
                resvalue = resvalueTwo;
                RandSwitch = -1;
            elseif resvalueTwo == resvalue % If both are the same then choose one at random
                RandSwitch = 2*round(rand())-1;         
            end

            % Compare to the current objective value
            if resvalue < objvalue % If an improvement is found:
                % Update objective value
                objvalue = resvalue; % Update objective value
                % Update Lower Design
                lowdes(:, ii) = RandSwitch*lowdes(:, ii);
                break % Get out of the most inner loop
            end
            % If no better design, then re-store lowdes
            lowdes(:, ii) = -1*lowdes(:, ii);
            lowdes(:, [ii jj]) = lowdes(:, [jj ii]); % Swap two columns
            
        end

    end

    % Evaluate if it is necessary to pass by the columns again
    thetadiff = abs(thetazero - objvalue);
    if thetadiff < 0.001 
        bflag = 2; % If there is no change in the solution then terminate
    else
        thetazero = objvalue; % If there is a change in the solution then continue and pass by the columns again
    end

end


end