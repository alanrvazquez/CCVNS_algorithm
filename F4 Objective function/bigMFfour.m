function [ fDvalue, CFV ] = bigMFfour(design, N, vecbigM, Selmintwo, mchoosetwo)

% Compute the f(D) values of a design. For more information, see
% Supplementary Section A.2.

% INPUTS:
% design        An N-by-m strength-3 design with coded levels -1 and +1 
%               (matrix).
% N             The run size of the design.
% vecbigM       The 1-by-nelF4 vector with the weight for each component of 
%               the F4 vector, where nelF4 is the number of elements in the 
%               F4 vector of a strength-3 design.
% mchoosetwo    The number of combinations of m in 2.
% Selmintwo     The mchoosetwo-by-m matrix of 2-factor subsets among m factors.
%
% OUTPUTS:
% fDvalue       The value of the f(D) function.
% CFV           The confounding frequency vector of length 4.
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

% Set up inputs and variables.---------------------------------------------
Int = TwoFIMat( design, Selmintwo, mchoosetwo, N );
nelF4 = ceil(N/16); % Number of elements in the F4 vector.
CFV = zeros(1, nelF4);
J4 = abs(Int'*Int); % Compute J4-characteristics.
J4 = J4 - diag(N*ones(1, size(J4, 1))); % Remove diagonal.

% Compute confounding frequency vector of a strength-3.-------------------- 
for clind=1:nelF4
    CFV(clind) = sum(sum(J4==N-(clind-1)*16))/6;
end

fDvalue = vecbigM*CFV'; % Objective value.---------------------------------
end

