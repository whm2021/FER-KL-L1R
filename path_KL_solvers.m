%%%%%%%%%%%%%%%%%%%%% KL divergence : running solvers %%%%%%%%%%%%%%%%%%%%%
% This script is called by the main.m script.
% It requires all simulation parameters to stored on 'param' variable.
% Variables A, y, lambda and x0 also need to be set.

%% Coordinate Descent - Hsieh2011
fprintf('CoD solver KL...\n')
[x_opt.CoD{i}, x_it.CoD{i}, primal_obj.CoD{i}, theta_opt.CoD{i}, stop_crit_it.CoD{i}, Iter.CoD{i}, time_it.CoD{i}] = CoD_KL_l1(A,y,lambda,para);
%% CoD + Screening
%% 情况：动态
fprintf('\n CoD solver KL + GAPScreening...\n')
precalc_GAP= KL_GAP_precalc(y, lambda, precalc);
[x_opt.CoD_GAP{i}, x_it.CoD_GAP{i}, primal_obj.CoD_GAP{i}, theta_opt.CoD_GAP{i}, stop_crit_it.CoD_GAP{i}, Iter.CoD_GAP{i}, time_it.CoD_GAP{i}, time_it_screen.CoD_GAP{i}, num_screen.CoD_GAP{i}] ...
    = CoD_KL_l1_GAP(A, y, lambda, para, precalc_GAP);

fprintf('\n CoD solver KL + G-GAPScreening...\n')
precalc_G_GAP= KL_G_GAP_precalc(A, y, lambda, precalc);
[x_opt.CoD_G_GAP{i}, x_it.CoD_G_GAP{i}, primal_obj.CoD_G_GAP{i}, theta_opt.CoD_G_GAP{i}, stop_crit_it.CoD_G_GAP{i}, Iter.CoD_G_GAP{i}, time_it.CoD_G_GAP{i}, time_it_screen.CoD_G_GAP{i}, num_screen.CoD_G_GAP{i}] ...
    = CoD_KL_l1_GAP(A, y, lambda, para, precalc_G_GAP);

fprintf('\n CoD solver KL + R-GAPScreening...\n')
precalc_G_GAP= KL_G_GAP_precalc(A, y, lambda, precalc);
[x_opt.CoD_R_GAP{i}, x_it.CoD_R_GAP{i}, primal_obj.CoD_R_GAP{i}, theta_opt.CoD_R_GAP{i}, stop_crit_it.CoD_R_GAP{i}, Iter.CoD_R_GAP{i}, time_it.CoD_R_GAP{i}, time_it_screen.CoD_R_GAP{i}, num_screen.CoD_R_GAP{i}] ...
    = CoD_KL_l1_R_GAP(A, y, lambda, para, precalc_G_GAP);

%% 情况：静态
fprintf('\n CoD solver KL + StaticScreening...\n')
precalc_GAP= KL_GAP_precalc(y, lambda, precalc);
if i==1
   [x_opt.CoD_Sta{i}, x_it.CoD_Sta{i}, primal_obj.CoD_Sta{i}, theta_opt.CoD_Sta{i}, stop_crit_it.CoD_Sta{i}, Iter.CoD_Sta{i}, time_it.CoD_Sta{i}]...
       = CoD_KL_l1(A,y,lambda,para);
else
   [x_opt.CoD_Sta{i}, x_it.CoD_Sta{i}, primal_obj.CoD_Sta{i}, theta_opt.CoD_Sta{i}, stop_crit_it.CoD_Sta{i}, Iter.CoD_Sta{i}, time_it.CoD_Sta{i}, time_it_screen.CoD_Sta{i}, num_screen.CoD_Sta{i}] ...
   = CoD_KL_l1_Sta(A, y, lambda, lambdas(i-1), x_opt.CoD_Sta{1, i-1}, theta_opt.CoD_Sta{1,i-1}, para, precalc_GAP);
end

fprintf('\n CoD solver KL + G-StaticScreening...\n')
precalc_G_GAP= KL_G_GAP_precalc(A, y, lambda, precalc);
if i==1
    [x_opt.CoD_G_Sta{i}, x_it.CoD_G_Sta{i}, primal_obj.CoD_G_Sta{i}, theta_opt.CoD_G_Sta{i}, stop_crit_it.CoD_G_Sta{i}, Iter.CoD_G_Sta{i}, time_it.CoD_G_Sta{i}]...
        = CoD_KL_l1(A,y,lambda,para);
else
    [x_opt.CoD_G_Sta{i}, x_it.CoD_G_Sta{i}, primal_obj.CoD_G_Sta{i}, theta_opt.CoD_G_Sta{i}, stop_crit_it.CoD_G_Sta{i}, Iter.CoD_G_Sta{i}, time_it.CoD_G_Sta{i}, time_it_screen.CoD_G_Sta{i}, num_screen.CoD_G_Sta{i}] ...
    = CoD_KL_l1_Sta(A, y, lambda, lambdas(i-1), x_opt.CoD_G_Sta{1, i-1}, theta_opt.CoD_G_Sta{1,i-1}, para, precalc_G_GAP);
end

%% 情况：混合筛选准则（静态+动态）
fprintf('\n CoD solver KL + Sta-GAPScreening...\n')
precalc_GAP= KL_GAP_precalc(y, lambda, precalc);
if i==1
   [x_opt.CoD_Sta_GAP{i}, x_it.CoD_Sta_GAP{i}, primal_obj.CoD_Sta_GAP{i}, theta_opt.CoD_Sta_GAP{i}, stop_crit_it.CoD_Sta_GAP{i}, Iter.CoD_Sta_GAP{i}, time_it.CoD_Sta_GAP{i}, time_it_screen.CoD_Sta_GAP{i}, num_screen.CoD_Sta_GAP{i}] ...
       = CoD_KL_l1_GAP(A, y, lambda, para, precalc_GAP);
else
   [x_opt.CoD_Sta_GAP{i}, x_it.CoD_Sta_GAP{i}, primal_obj.CoD_Sta_GAP{i}, theta_opt.CoD_Sta_GAP{i}, stop_crit_it.CoD_Sta_GAP{i}, Iter.CoD_Sta_GAP{i}, time_it.CoD_Sta_GAP{i}, time_it_screen.CoD_Sta_GAP{i}, num_screen.CoD_Sta_GAP{i}] ...
       = CoD_KL_l1_Sta_GAP(A, y, lambda, lambdas(i-1), x_opt.CoD_Sta_GAP{1, i-1}, theta_opt.CoD_Sta_GAP{1,i-1}, para, precalc_GAP);
end

fprintf('\n CoD solver KL + G-Sta-GAPScreening...\n')
precalc_G_GAP= KL_G_GAP_precalc(A, y, lambda, precalc);
if i==1
    [x_opt.CoD_G_Sta_GAP{i}, x_it.CoD_G_Sta_GAP{i}, primal_obj.CoD_G_Sta_GAP{i}, theta_opt.CoD_G_Sta_GAP{i}, stop_crit_it.CoD_G_Sta_GAP{i}, Iter.CoD_G_Sta_GAP{i}, time_it.CoD_G_Sta_GAP{i}, time_it_screen.CoD_G_Sta_GAP{i}, num_screen.CoD_G_Sta_GAP{i}] ...
        = CoD_KL_l1_GAP(A, y, lambda, para, precalc_G_GAP);
else
    [x_opt.CoD_G_Sta_GAP{i}, x_it.CoD_G_Sta_GAP{i}, primal_obj.CoD_G_Sta_GAP{i}, theta_opt.CoD_G_Sta_GAP{i}, stop_crit_it.CoD_G_Sta_GAP{i}, Iter.CoD_G_Sta_GAP{i}, time_it.CoD_G_Sta_GAP{i}, time_it_screen.CoD_G_Sta_GAP{i}, num_screen.CoD_G_Sta_GAP{i}] ...
        = CoD_KL_l1_Sta_GAP(A, y, lambda, lambdas(i-1), x_opt.CoD_G_Sta_GAP{1, i-1}, theta_opt.CoD_G_Sta_GAP{1,i-1}, para, precalc_G_GAP);
end

















