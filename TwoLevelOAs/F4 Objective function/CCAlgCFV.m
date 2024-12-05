function [objvalue, lowdes] = CCAlgCFV(lowdes, m, UTU, UpInt, Selmintwo, nelF4, vecbigM, Ndiag, PossibJ4s, mchoosetwo, N) 

% Column Change (CC) algorithm for the optimizing the f(D) function based on 
% the F4 vector. 
%
% INPUTS:
% lowdes        An N-by-m upper design with coded levels -1 and +1 (matrix).
% m             The number of factors.
% mchoosetwo    The number of combinations of m in 2.
% UpInt         The N-by-mchoosetwo two-factor (2FI) interaction matrix 
%               for the upper design.
% UTU           The mchoosetwo-by-mchoosetwo matrix containing the product 
%               UpInt'*UpInt.
% Selmintwo     The mchoosetwo-by-m matrix of 2-factor subsets among m factors.
% nelF4         The number of elements in the F4 vector of a strength-3 design.
% vecbigM       The 1-by-nelF4 vector with the weight for each component of the
%               F4 vector.
% Ndiag         The 2N-by-2N identity matrix.
% PossibJ4s     Possible values of the J4-characteristics.
% N             The run size of the upper or lower design.
%
% OUTPUTS:
% objvalue      The best objective value found.
% lowdes        An N-by-m matrix containing the improved plan for the lower 
%               design.
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

%====SET UP CURRENT OBJECTIVE VALUE========================================
LowInt = TwoFIMat(lowdes, Selmintwo, mchoosetwo, N); % 2FI matrix.
BB = LowInt'*LowInt;
J4 = abs(UTU + BB - Ndiag); % Compute J4-characteristics using matrix 
                            % computations.
J4vec = J4(:);
CFV = histc(J4vec, PossibJ4s)'/6;
CFV = CFV(nelF4:-1:1); % Compute Confounding Frequency (CF) vector.
objvalue = vecbigM*(CFV)'; % Compute objective value.

%====FLAGS TO STOP THE ALGORITHM===========================================
thetazero=0;
bflag=0;

%====START CC ALGORITHM====================================================
while bflag < 1
% Step 2 Improve
    for ii = 1:m

        % Sign switch column 'ii'. Compute updating formula.---------------
        [CFVCont, SelectFact] = SignSwitchImpact(ii, LowInt, UpInt, Selmintwo, nelF4, PossibJ4s );
        % Remove J4-characteristics of column 'ii' and add the 
        % J4-characteristics of column '-ii'.
        resvalue = vecbigM*(CFV - CFVCont)'; 

        if resvalue < objvalue % If improvement.
            objvalue = resvalue; % Update best objective value.
            lowdes(:, ii) = -1*lowdes(:, ii); % Update best lower design.
            CFV = CFV - CFVCont; % Update CF vector.
            LowInt(:, SelectFact) = -1*LowInt(:, SelectFact); % Update 2FI
                                                              % matrix for 
                                                              % the lower
                                                              % design
            continue % Do not execute the following and go to 'ii+1'.
            
        end

        % Swap column 'ii' and the columns to its right.-------------------
        st1 = ii+1;
        for jj = st1:m
   
            Auxlowdes = lowdes;
            % First. Compute the contributions of each move.---------------
            % Swapping columns 'ii' and 'jj' or columns '-ii' and 'jj'.
            
            % Compute the CF vector for a concatenated design constructed 
            % from a lower design formed by by swapping columns 'ii' and 'jj'.
            % Report also the 2FI matrix for this plan of the lower design.
            % AuxLowInt: plan for the lower design formed by swapping columns 
            % 'ii' and 'jj'.
            [ CFVAux, AuxLowInt ] = ChangeColImpact(ii, jj, Auxlowdes, LowInt, UTU, Selmintwo, nelF4, Ndiag, PossibJ4s);
            
            % Compute the CF vector for a concatenated design constructed 
            % from a lower design formed by by swapping columns '-ii' and 'jj'.
            % Report also the 2FI matrix for this plan of the lower design.
            % SelectFact: A 1-by-mchoosetwo logical vector that determines 
            % the positions of the elements involving factor 'ii' in the lower 
            % design, where mchoosetwo is the number of combinations of m 
            % factors in two.
            [ CFVNeg, SelectFact ] = SignSwitchImpact(ii, AuxLowInt, UpInt, Selmintwo, nelF4, PossibJ4s );
            
            resvalue =  vecbigM*CFVAux'; % Objective value when swapping 
                                         % columns 'ii' and 'jj'.
            resvalueTwo = vecbigM*(CFVAux - CFVNeg)'; % Objective value 
                                                      % when swapping columns 
                                                      % '-ii' and 'jj'.
            % Swap columns in the design.----------------------------------
            Auxlowdes(:, [ii jj]) = Auxlowdes(:, [jj ii]); 
            
            % Determine which one is the best change.----------------------
            RandSwitch = 1;
            if resvalueTwo < resvalue
                 resvalue = resvalueTwo; 
                 RandSwitch = -1;
            elseif resvalueTwo == resvalue % If both are the same then choose 
                                           % one at random.
                 RandSwitch = 2*round(rand())-1; 
            end

            % Compare to the current objective value.----------------------
            if resvalue < objvalue % If improvement.
                objvalue = resvalue; % Update best objective value.
    
                Auxlowdes(:, ii) = RandSwitch*Auxlowdes(:, ii); 
                lowdes = Auxlowdes; % Update best lower design.
                
                AuxLowInt(:, SelectFact) = RandSwitch*AuxLowInt(:, SelectFact);
                LowInt = AuxLowInt; % Update 2FI matrix of the lower design.
                
                % Update CF vector.
                CFV = CFVAux + ( (RandSwitch - 1)/2 )*CFVNeg;
                break % Get out of the most inner loop.                 
            end

        end % end inner for.

    end % end outer for.

    % Evaluate if it is necessary to pass by the columns again-------------
    thetadiff = abs(thetazero - objvalue);
    if thetadiff < 0.001 
        bflag = 2; % If there is no change in the solution then terminate.
    else
        thetazero = objvalue; % If there is a change in the solution then 
                              % continue and pass by the columns again.
    end

end % end while.


end