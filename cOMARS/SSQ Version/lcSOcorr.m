function [ res ] = lcSOcorr( design, ncombtwo, nm, N)

% lcSOcorr Sum of square correlations between pairs of second-order
%           effect columns.

Intmat = TwoFIMat( design, ncombtwo, nm, N ); % Two-factor interaction matrix
correlations_int = corr(Intmat);
triu_corr_int = triu(correlations_int, 1);
res = triu_corr_int(:)'*triu_corr_int(:);

end