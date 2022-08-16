clear;clear all
% ===================== 将所有用到的文件夹都添加到路径 =====================
addpath('./solvers/')
addpath('./screening_rules/')
addpath('./datasets/')
load('test_data.mat')
% load('data_20newsgroup.mat')
% load('data_AllBooks.mat')
% load('data_Encyclopedia.mat')
% load('data_MNIST.mat')
% ============================= 固定随机数 ================================
rng_seed = 10; % 0 for no seed （设置随机变量的种子，seed一旦确定后，每次产生的随机值都是一样的了）
if rng_seed
    rng(rng_seed) % 固定随机数必须使用该命令
    fprintf('\n\n /!\\/!\\ RANDOM SEED ACTIVATED /!\\/!\\\n\n'); 
end
% ============================= 预先计算值 ================================
Ap = pinv(full(A));
precalc.NormAp1 = sum(abs(Ap))';
precalc.NormA = sqrt(sum(A(y~=0,:)'.^2,2));
precalc.NormA1 = norm(A, 1);
precalc.NormA11 = sum(abs(A));
% ============================= 参数设置 ==================================
para.tol = 1e-7; %对偶间隙<tol
para.max_iter = 1e5; %最大迭代次数
para.epsilon = 1e-6; %避免KL散度中奇异值为0的情况
para.tol_r = 1e-8;
para.tol_gap = 0.05; %设置GAP筛选到一定程度就终止了
% 迭代参数lambda设置
t = 1:1:100;
T = 100;
deta = 4;
inter = 10.^-((deta*t)./(T-1));
lambda_max = max(A.'*(y-para.epsilon))/para.epsilon;
lambdas = lambda_max.*inter;
% ============================= 存储变量设置 ===============================
% 原问题最优解
x_opt.CoD = cell(1, T);
x_opt.CoD_GAP = cell(1, T);
x_opt.CoD_G_GAP = cell(1, T);
x_opt.CoD_R_GAP = cell(1, T);
x_opt.CoD_Sta = cell(1, T);
x_opt.CoD_G_Sta = cell(1, T);
x_opt.CoD_Sta_GAP = cell(1, T);
x_opt.CoD_G_Sta_GAP = cell(1, T);
% 每次迭代得到的可行解
x_it.CoD = cell(1, T);
x_it.CoD_GAP = cell(1, T);
x_it.CoD_G_GAP = cell(1, T);
x_it.CoD_R_GAP = cell(1, T);
x_it.CoD_Sta = cell(1, T);
x_it.CoD_G_Sta = cell(1, T);
x_it.CoD_Sta_GAP = cell(1, T);
x_it.CoD_G_Sta_GAP = cell(1, T);
% 原问题目标函数值
primal_obj.CoD = cell(1, T);
primal_obj.CoD_GAP = cell(1, T);
primal_obj.CoD_G_GAP = cell(1, T);
primal_obj.CoD_R_GAP = cell(1, T);
primal_obj.CoD_Sta = cell(1, T);
primal_obj.CoD_G_Sta = cell(1, T);
primal_obj.CoD_Sta_GAP = cell(1, T);
primal_obj.CoD_G_Sta_GAP = cell(1, T);
% 对偶问题最优解
theta_opt.CoD = cell(1, T);
theta_opt.CoD_GAP = cell(1, T);
theta_opt.CoD_G_GAP = cell(1, T);
theta_opt.CoD_R_GAP = cell(1, T);
theta_opt.CoD_Sta = cell(1, T);
theta_opt.CoD_G_Sta = cell(1, T);
theta_opt.CoD_Sta_GAP = cell(1, T);
theta_opt.CoD_G_Sta_GAP = cell(1, T);
% 终止准则(gap的变化规律)
stop_crit_it.CoD = cell(1, T);
stop_crit_it.CoD_GAP = cell(1, T);
stop_crit_it.CoD_G_GAP = cell(1, T);
stop_crit_it.CoD_R_GAP = cell(1, T);
stop_crit_it.CoD_Sta = cell(1, T);
stop_crit_it.CoD_G_Sta = cell(1, T);
stop_crit_it.CoD_Sta_GAP = cell(1, T);
stop_crit_it.CoD_G_Sta_GAP = cell(1, T);
% 迭代次数
Iter.CoD = cell(1, T);
Iter.CoD_GAP = cell(1, T);
Iter.CoD_G_GAP = cell(1, T);
Iter.CoD_R_GAP = cell(1, T);
Iter.CoD_Sta = cell(1, T);
Iter.CoD_G_Sta = cell(1, T);
Iter.CoD_Sta_GAP = cell(1, T);
Iter.CoD_G_Sta_GAP = cell(1, T);
% 每次迭代的时间
time_it.CoD = cell(1, T);
time_it.CoD_GAP = cell(1, T);
time_it.CoD_G_GAP = cell(1, T);
time_it.CoD_R_GAP = cell(1, T);
time_it.CoD_Sta = cell(1, T);
time_it.CoD_G_Sta = cell(1, T);
time_it.CoD_Sta_GAP = cell(1, T);
time_it.CoD_G_Sta_GAP = cell(1, T);
% 筛选时间
time_it_screen.CoD_GAP = cell(1,T);
time_it_screen.CoD_G_GAP = cell(1,T);
time_it_screen.CoD_R_GAP = cell(1,T);
time_it_screen.CoD_Sta = cell(1,T);
time_it_screen.CoD_G_Sta = cell(1,T);
time_it_screen.CoD_Sta_GAP = cell(1,T);
time_it_screen.CoD_G_Sta_GAP = cell(1,T);
% 筛选数量
num_screen.CoD_GAP= cell(1, T);
num_screen.CoD_G_GAP= cell(1, T);
num_screen.CoD_R_GAP= cell(1, T);
num_screen.CoD_Sta= cell(1, T);
num_screen.CoD_G_Sta= cell(1, T);
num_screen.CoD_Sta_GAP= cell(1, T);
num_screen.CoD_G_Sta_GAP= cell(1, T);
%% =========================main loop======================================
for i = 1 : length(lambdas)
    lambda = lambdas(i);
    fprintf('\n ---- Regularization parameter %d / %d ----\n',i, length(lambdas))
    path_KL_solvers
end
save('result_test', 'x_opt', 'x_it', 'primal_obj', 'theta_opt', 'stop_crit_it','Iter','time_it','time_it_screen','num_screen')

%% 计算最终的结果
Compute_final_result













