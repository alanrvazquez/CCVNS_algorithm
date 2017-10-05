function [ updateCFV, SelectFact ] = SignSwitchImpact(ii, LowInt, UpInt, Selmintwo, nelF4, PossibJ4s )

% Calculate the J4-characteristics over 4-factor sets involving column 'ii'.
% For more information, see Supplementary Section A.2.
%
% INPUTS:
% ii            The column to be changed (scalar).
% LowInt        The N-by-mchoosetwo two-factor (2FI) interaction matrix for 
%               the lower design, where mchoosetwo is the number of 
%               combinations of m factors in 2.
% UpInt         The N-by-mchoosetwo two-factor (2FI) interaction matrix for 
%               the upper design.
% Selmintwo     The mchoosetwo-by-m matrix of 2-factor subsets among m factors.
% nelF4         The number of elements in the F4 vector of a strength-3 design.
% PossibJ4s     Possible values of the J4-characteristics.
%
% OUTPUTS:
% updateCFV     Confounding frequency (CF) vector computed as the
%               difference between the CF vector of a concatenated design
%               constructed from a plan of the lower design that includes 
%               column 'ii' and the CF vector of a concatenated design 
%               constructed from a plan of the lower design that includes
%               column '-ii'. The CF vectors are computed only for the 
%               J4-characteristics that involve column 'ii'.
% SelectFact    A 1-by-mchoosetwo logical vector that determines the positions 
%               of the elements involving factor 'ii' in the lower design,
%               where mchoosetwo is the number of combinations of m factors
%               in two.
%
% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management
%==========================================================================

% Logical vectors to select the elements involving factor ii in the 2FI 
% matrices.----------------------------------------------------------------
SelectFact = or(Selmintwo(:, 1) == ii, Selmintwo(:, 2) == ii);
NotSelFact = ~SelectFact;

% Compute matrix El'*Gl. For details, see equation 3 in Supplementary
% Section A.2.-------------------------------------------------------------
IntFacti = LowInt(:, SelectFact')'*LowInt(:, NotSelFact'); 
% Compute matrix Eu'*Gu.
UTU = UpInt(:, SelectFact')'*UpInt(:, NotSelFact'); 

% Compute J4-characteristics.----------------------------------------------
J4Pos = abs(UTU + IntFacti); % NO sign switch.
J4Neg = abs(UTU - IntFacti); % Sign switch.
J4Posvec = J4Pos(:);
J4Negvec = J4Neg(:);

% Compute the CF vector of length 4 over J4-characteristics involving column 
% '-ii'.-------------------------------------------------------------------
% Remove the J4-characteristics of column 'ii' and add the J4-characteristics 
% of column '-ii'.
CFVCont = histc(J4Posvec, PossibJ4s)' - histc(J4Negvec, PossibJ4s)';
CFVCont = CFVCont(nelF4:-1:1);
updateCFV = CFVCont/3;
end

