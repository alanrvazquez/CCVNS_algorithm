function [ CFVCalc , LowInt ] = ChangeColImpact(ii, jj, deslower, LowInt, UTU, Selmintwo, nelF4, Ndiag, PossibJ4s )

% Calculate the J4-characteristics over 4-factor sets involving column 'ii'.
% For more information, see Supplementary Section A.2.
%
% INPUTS:
% ii, jj        The columns to be swapped (scalars).
% deslower      The N-by-m lower design with coded levels +1 and -1.
% LowInt        The N-by-mchoosetwo two-factor (2FI) interaction matrix for 
%               the lower design, where mchoosetwo is the number of 
%               combinations of m factors in 2.
% UTU           The mchoosetwo-by-mchoosetwo matrix containing the product 
%               UpInt'*UpInt, where UpInt is the N-by-n two-factor (2FI) 
%               interaction matrix for the upper design.
% Selmintwo     The mchoosetwo-by-m matrix of 2-factor subsets among m
%               factors.
% nelF4         The number of elements in the F4 vector of a strength-3 design.
% Ndiag         The 2N-by-2N identity matrix.
% PossibJ4s     Possible values of the J4-characteristics.
%
% OUTPUTS:
% CFVCalc       Confounding frequency (CF) vector of a concatenated design
%               construted from a plan of the lower design that swaps
%               columns 'ii' and 'jj'.
% LowInt        The N-by-n two-factor (2FI) interaction matrix for the plan
%               of the lower design that swaps columns 'ii' and 'jj'. 
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

% Compute the J4-characteristics of the new plan of the lower design.------ 
zz = deslower(:, ii).*deslower(:, jj); 
SelectColumns = xor(sum(Selmintwo == ii, 2) ~=0, sum(Selmintwo == jj, 2) ~=0)';
% Fast update of the 2FI interaction of the lower design.
LowInt(:, SelectColumns) = bsxfun(@times,zz,LowInt(:, SelectColumns)); 
Matt = abs(UTU + LowInt'*LowInt) - Ndiag; % J4-characteristics.

% Compute the CF vector for the new plan of the lower design.--------------
MattVec = Matt(:);
CFVCalc = histc(MattVec, PossibJ4s)';
CFVCalc = CFVCalc(nelF4:-1:1);
CFVCalc = CFVCalc/6;
end

