function [ Bfour ] = Bfour( D, N, m)

% Bfour Compute B4 value of two-level designs using Butler's moment matrix. 
% For details see: 
% Butler, N. A. (2003). Minimum aberration construction results for
% nonregular two-level factorial designs. Biometrika, 90:891-898.
%
% INPUTS:
% D   The N-by-m design with coded levels -1 and +1 (matrix).
% N   The run size of the design.
% m   The number of factors.
%
% OUTPUTS:
% Bfour  B4 value of the design.
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

H = D*D'; % Compute Moment matrix------------------------------------------

% Fourth moment------------------------------------------------------------
HvectoF = H(:).^4;
Mfour = (N^(-2))*(sum(HvectoF)); 
Bfour = (Mfour - m*(3*m-2))/24; % Translate into from M4 into B4-----------

end

