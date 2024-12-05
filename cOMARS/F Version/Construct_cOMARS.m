%% =====CONSTRUCT A CONCATENATED DESIGN=====================================

% This script shows how to construct a cOMARS design that minimizes 
% the sum of squared correlations between pairs of second-order effect
% columns.

% AUTHOR: 
% Alan R. Vazquez
% Department of Industrial Engineering
% Tecnologico de Monterrey

% Set the number of factors.-----------------------------------------------
m = 19;
id_one = 1;
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
[obj_val, design] = CCVNSAlg(upper_design, lower_design, max_iterations);
%delete(pobj); % Shutdown parallel pool.

% Select best design.------------------------------------------------------
select_best = lab_vec(min(obj_val) == obj_val);
best_design = design{select_best(1)};

% Construct cOMARS design.-------------------------------------------------
cOMARS = [best_design; -1*best_design; zeros(nzero, m)];
disp(cOMARS); % Display cOMARS design.

%% Save cOMARS design.------------------------------------------------------
save_path = strcat("cOMARS/cOMARS_m_", strm, "_id_", strid_one, "_", strid_two, ".csv");
save_path_obj = strcat("cOMARS_m_", strm, "_id_", strid_one, "_", strid_two,  ".mat");

writematrix(cOMARS, save_path)
save(save_path_obj, "design", "obj_val")