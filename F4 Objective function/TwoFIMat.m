function [ Matt ] = TwoFIMat( design, Selmintwo, mchoosetwo, N )

% Compute the two-factor interaction matrix of a design.
% Note: This function is faster than the function to construct a two-factor 
% interaction matrix in Matlab, xfsx(M, 'i').
%
% INPUTS:
% design        An N-by-m design with coded levels -1 and +1 (matrix).
% mchoosetwo    The number of combinations of m in 2.
% Selmintwo     The mchoosetwo-by-m matrix of 2-factor subsets among m factors.
% N             The run size of the design.
%
% OUTPUTS:
% Matt          An N-by-mchoosetwo two-factor interaction matrix. 
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

Matt = zeros(N, mchoosetwo);
for ii = 1:mchoosetwo   
    Matt(:, ii) = design(:, Selmintwo(ii, 1)).*design(:, Selmintwo(ii, 2));
end

end

