function [objvalue, lowdes] = CCAlgBfour(desupper, lowdes, m)

% Column Change (CC) algorithm for optimizing the B4 value. 
%
% INPUTS:
% desupper  An N-by-m upper design with coded levels -1 and +1 (matrix).
% lowdes    An N-by-m lower design with coded levels -1 and +1 (matrix).
% m         The number of factors.
%
% OUTPUTS:
% objvalue  The best objective value found.
% lowdes    An N-by-m matrix containing the improved plan for the lower 
%           design.
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

%====SET UP CURRENT OBJECTIVE VALUE========================================
Hmat = desupper*lowdes';
mattAux = Hmat(:);
threeauxval = mattAux.^2;
objvalue = threeauxval'*threeauxval; 

%====FLAGS TO STOP THE ALGORITHM===========================================
thetazero=0;
bflag=0;

%====START CC ALGORITHM====================================================
while bflag < 1
    for ii = 1:m

        % Sign switch column 'ii'. Compute updating formula.---------------
        FoldHmat = Hmat - 2*desupper(:, ii)*lowdes(:, ii)'; 
        vecAux = FoldHmat(:).^2;
        resvalue = vecAux'*vecAux;

        if resvalue < objvalue % If improvement.
            objvalue = resvalue; % Update best objective value.
            lowdes(:, ii) = -1*lowdes(:, ii); % Update best lower design.
            Hmat = FoldHmat; % Update Hmat.
            continue % Do not execute the following and go to 'ii+1'.
        end

        % Swap column 'ii' and the columns to its right.-------------------
        st1 = ii+1;
        for jj = st1:m
   
            Auxlowdes = lowdes; 
            Auxlowdes(:, [ii jj]) = Auxlowdes(:, [jj ii]); % Swap two columns.
            
            % Compute the updating formulas.------------------------------- 
            MatValueOne = desupper*Auxlowdes'; % Just swap two columns 
                                               % (Positive version).
            MatValueTwo = MatValueOne - 2*desupper(:, ii)*Auxlowdes(:, ii)'; % Swap and change sings of column 'ii' 
                                                                             % (Negative version).
            
            % Evaluate positive and negative versions.---------------------
            VecValueOne = MatValueOne(:).^2; 
            VecValueTwo =  MatValueTwo(:).^2; 
            resvalue = VecValueOne'*VecValueOne; 
            resvalueTwo = VecValueTwo'*VecValueTwo; 

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
                Hmat = desupper*lowdes'; % Update Hmat.
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