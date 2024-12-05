%% =====CONSTRUCT A CONCATENATED DESIGN=====================================

% This script shows how to construct a cOMARS design that minimizes 
% the sum of squared correlations between pairs of second-order effect
% columns.

% AUTHOR: 
% Alan R. Vazquez
% Department of Industrial Engineering
% Tecnologico de Monterrey

% Set the number of factors.-----------------------------------------------
m = 20;
id_one = 2;
id_two = 2;

% Set number of zeros in the cOMARS design.--------------------------------
nzero = 1; 

% Set maximum number of iterations.----------------------------------------
max_iterations = 10;

% ==== LOAD CONFERENCE MATRICES ===========================================

% Select a conference design.
strm = num2str(m);
strid_one = num2str(id_one);
sel_des_one = strcat("m", strm, "_", strid_one);

strm = num2str(m);
strid_two = num2str(id_two);
sel_des_two = strcat("m", strm, "_", strid_two);

% Set the parent designs.--------------------------------------------------
upper_design = readmatrix("Conference Designs.xlsx", Sheet= sel_des_one); 
lower_design = readmatrix("Conference Designs.xlsx", Sheet= sel_des_two); 

% Set other parameters.----------------------------------------------------
lab_vec = 1:max_iterations;
rng(81479058); % Set seed.

% Parallel version of the algorithm.---------------------------------------
%pobj = parpool('local'); % Create parallel pool.
[SSQcor_val, design] = CCVNSAlg(upper_design, lower_design, max_iterations);
%delete(pobj); % Shutdown parallel pool.

% Select best design.------------------------------------------------------
select_best = lab_vec(min(SSQcor_val) == SSQcor_val);
best_design = design{select_best(1)};

% Construct cOMARS design.-------------------------------------------------
cOMARS = [best_design; -1*best_design; zeros(nzero, m)];

None = nchoosek(m, 2);
NSelTwo = nchoosek(1:m, 2);
ssq_comars = lcSOcorr( cOMARS, NSelTwo, None, size(cOMARS,1));

disp(cOMARS); % Display cOMARS design.
disp(ssq_comars); % Display sum of squared correlations

% Save cOMARS design.------------------------------------------------------
save_path = strcat("cOMARS/cOMARS_m_", strm, "_id_2.csv");
writematrix(cOMARS, save_path)

save_path_obj = strcat("cOMARS_m_", strm, "_id_2.mat");
save(save_path_obj, "design", "SSQcor_val")