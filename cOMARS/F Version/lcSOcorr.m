function [ res ] = lcSOcorr( design, ncombtwo, nm, N, vecbigM, nelLambda, Lambda)

% lcSOcorr Sum of square correlations between pairs of second-order
%           effect columns.

Intmat = TwoFIMat( design, ncombtwo, nm, N ); % Two-factor interaction matrix
CFV = zeros(1, nelLambda);
J4 = abs(Intmat'*Intmat); % Compute J4-characteristics.
J4 = J4 - diag(N*ones(1, size(J4, 1))); % Remove diagonal.

% Compute confounding frequency vector of a strength-3.-------------------- 
for clind=1:nelLambda
    CFV(clind) = sum(sum(J4==Lambda(clind)));
end

res = vecbigM*CFV';

end