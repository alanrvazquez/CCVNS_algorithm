%% =====CONSTRUCT A CONCATENATED DESIGN=====================================

% This script shows how to construct an concatenated designs that optimizes 
% the B4 value.

% AUTHOR: 
% Alan Vazquez-Alcocer
% University of Antwerp
% Department of Engineering Management

%% ====Example: Construct a concatenated design with 64 runs and 17 factors

% Load strength-3 orthogonal arrays with 32 runs and 16 factors.
% There are five arrays in total. These are obtained from the complete
% enumeration of Schoen, E. D., Eendebak P. T., and Nguyen M. V. M. (2010).
% Complete enumeration of pure-level and mixed-level orthogonal arrays.
% Journal of Combinatorial Designs, 18:-123-140.

% For more details visit: http://www.pietereendebak.nl/oapackage/ 
load('OAN32K16') % A 32-by-16-by-5 array containing five orthogonal arrays 
                 % with coded levels 0 and 1.

% Select two parent designs.-----------------------------------------------
% Code designs to -1s and 1s.
upper_design = 2*OAN32K16(:, :, 4)-1; % Design 4 in the catalog.
lower_design = 2*OAN32K16(:, :, 5)-1; % Design 5 in the catalog.

% Set maximum number of iterations.----------------------------------------
max_iterations = 40;

% Set other parameters.----------------------------------------------------
lab_vec = 1:max_iterations;
rng(442); % Set seed.
%% ====APPLY CC/VNS ALGORITHM TO OPTIMIZE CONCATENATED DESIGN==============
% Nowadays, almost all standard computers have multi-core processors. If
% this is the case for your computer, you might want to use the parallel
% version of the algorithm. It is more efficient as it distributes the
% interations into your computer processors. 

% Standard version of the algorithm---------------------------------------
tic
[Bfour_val, design] = CCVNSBfour(upper_design, lower_design, max_iterations );
toc
% Select best design.------------------------------------------------------
select_best = lab_vec(min(Bfour_val) == Bfour_val);
best_design = design{select_best(1)};
disp(best_design); % Print design. 
disp(Bfour_val(select_best(1))); % B4 value of the best design found over all 
                              % iterations.
% Compute the F4 vector and the B4 value of the design.--------------------
result = F4(best_design);
disp(result{1}); % F4 vector.
disp(result{2}); % B4 value.                              
% If needed, the complete list of designs can be evaluated according to
% other criteria so that the best (admissible) design can be chosen.

%% Parallel version of the algorithm.--------------------------------------
pobj = parpool('local'); % Create parallel pool.
[Bfour_val, design] = Parallel_CCVNSBfour(upper_design, lower_design, max_iterations );
delete(pobj); % Shutdown parallel pool.

% Select best design.------------------------------------------------------
select_best = lab_vec(min(Bfour_val) == Bfour_val);
best_design = design{select_best(1)};
disp(best_design); % Print design.  
disp(Bfour_val(select_best(1))); % B4 value of the best design found over all 
                              % iterations.
% Compute the F4 vector and the B4 value of the design.--------------------
result = F4(best_design);
disp(result{1}); % F4 vector.
disp(result{2}); % B4 value.                              
% If needed, the complete list of designs can be evaluated according to
% other criteria so that the best (admissible) design can be chosen.



