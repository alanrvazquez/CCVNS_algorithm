function [ Matt ] = TwoFIMat( design, ncombtwo, nm, N )
%TwoFIMat Compute the two-factor interaction model matrix
%   This is a faster alternative to the 'x2fx' function

% design: design under study coded to -1 and 1
% ncomtwo: matrix containing the combinations of m, number of factors, and
% 2
% nm: number of combinations of n in 2
% N number of runs of design

Matt = zeros(N, nm);
for ii = 1:nm
    Matt(:, ii) = design(:, ncombtwo(ii, 1)).*design(:, ncombtwo(ii, 2));  
end

end

