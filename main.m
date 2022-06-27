clear;clear all
% ===================== �������õ����ļ��ж���ӵ�·�� =====================
addpath('./solvers/')
addpath('./screening_rules/')
addpath('./datasets/')
%load('test_data.mat')
% load('data_20newsgroup.mat')
load('data_AllBooks.mat')
% load('data_Encyclopedia.mat')
%load('data_MNIST.mat')
% ============================= �̶������ ================================
rng_seed = 10; % 0 for no seed ������������������ӣ�seedһ��ȷ����ÿ�β��������ֵ����һ�����ˣ�
if rng_seed
    rng(rng_seed) % �̶����������ʹ�ø�����
    fprintf('\n\n /!\\/!\\ RANDOM SEED ACTIVATED /!\\/!\\\n\n'); 
end
% ============================= Ԥ�ȼ���ֵ ================================
Ap = pinv(full(A));
precalc.NormAp1 = sum(abs(Ap))';
precalc.NormA = sqrt(sum(A(y~=0,:)'.^2,2));
precalc.NormA1 = norm(A, 1);
precalc.NormA11 = sum(abs(A));
% ============================= �������� ==================================
para.tol = 1e-7; %��ż��϶<tol
para.max_iter = 1e5; %����������
para.epsilon = 1e-6; %����KLɢ��������ֵΪ0�����
para.tol_r = 1e-8;
para.tol_gap = 0.05; %����GAPɸѡ��һ���̶Ⱦ���ֹ��
% ��������lambda����
t = 1:1:100;
T = 100;
deta = 4;
inter = 10.^-((deta*t)./(T-1));
lambda_max = max(A.'*(y-para.epsilon))/para.epsilon;
lambdas = lambda_max.*inter;
% ============================= �洢�������� ===============================
% ԭ�������Ž�
x_opt.CoD = cell(1, T);
x_opt.CoD_GAP = cell(1, T);
x_opt.CoD_G_GAP = cell(1, T);
x_opt.CoD_R_GAP = cell(1, T);
x_opt.CoD_Dome_GAP = cell(1, T);
x_opt.CoD_Sta = cell(1, T);
x_opt.CoD_G_Sta = cell(1, T);
x_opt.CoD_Dome_Sta = cell(1, T);
x_opt.CoD_Sta_GAP = cell(1, T);
x_opt.CoD_G_Sta_GAP = cell(1, T);
% ÿ�ε����õ��Ŀ��н�
x_it.CoD = cell(1, T);
x_it.CoD_GAP = cell(1, T);
x_it.CoD_G_GAP = cell(1, T);
x_it.CoD_R_GAP = cell(1, T);
x_it.CoD_Dome_GAP = cell(1, T);
x_it.CoD_Sta = cell(1, T);
x_it.CoD_G_Sta = cell(1, T);
x_it.CoD_Dome_Sta = cell(1, T);
x_it.CoD_Sta_GAP = cell(1, T);
x_it.CoD_G_Sta_GAP = cell(1, T);
% ԭ����Ŀ�꺯��ֵ
primal_obj.CoD = cell(1, T);
primal_obj.CoD_GAP = cell(1, T);
primal_obj.CoD_G_GAP = cell(1, T);
primal_obj.CoD_R_GAP = cell(1, T);
primal_obj.CoD_Dome_GAP = cell(1, T);
primal_obj.CoD_Sta = cell(1, T);
primal_obj.CoD_G_Sta = cell(1, T);
primal_obj.CoD_Dome_Sta = cell(1, T);
primal_obj.CoD_Sta_GAP = cell(1, T);
primal_obj.CoD_G_Sta_GAP = cell(1, T);
% ��ż�������Ž�
theta_opt.CoD = cell(1, T);
theta_opt.CoD_GAP = cell(1, T);
theta_opt.CoD_G_GAP = cell(1, T);
theta_opt.CoD_R_GAP = cell(1, T);
theta_opt.CoD_Dome_GAP = cell(1, T);
theta_opt.CoD_Sta = cell(1, T);
theta_opt.CoD_G_Sta = cell(1, T);
theta_opt.CoD_Dome_Sta = cell(1, T);
theta_opt.CoD_Sta_GAP = cell(1, T);
theta_opt.CoD_G_Sta_GAP = cell(1, T);
% ��ֹ׼��(gap�ı仯����)
stop_crit_it.CoD = cell(1, T);
stop_crit_it.CoD_GAP = cell(1, T);
stop_crit_it.CoD_G_GAP = cell(1, T);
stop_crit_it.CoD_R_GAP = cell(1, T);
stop_crit_it.CoD_Dome_GAP = cell(1, T);
stop_crit_it.CoD_Sta = cell(1, T);
stop_crit_it.CoD_G_Sta = cell(1, T);
stop_crit_it.CoD_Dome_Sta = cell(1, T);
stop_crit_it.CoD_Sta_GAP = cell(1, T);
stop_crit_it.CoD_G_Sta_GAP = cell(1, T);
% ��������
Iter.CoD = cell(1, T);
Iter.CoD_GAP = cell(1, T);
Iter.CoD_G_GAP = cell(1, T);
Iter.CoD_R_GAP = cell(1, T);
Iter.CoD_Dome_GAP = cell(1, T);
Iter.CoD_Sta = cell(1, T);
Iter.CoD_G_Sta = cell(1, T);
Iter.CoD_Dome_Sta = cell(1, T);
Iter.CoD_Sta_GAP = cell(1, T);
Iter.CoD_G_Sta_GAP = cell(1, T);
% ÿ�ε�����ʱ��
time_it.CoD = cell(1, T);
time_it.CoD_GAP = cell(1, T);
time_it.CoD_G_GAP = cell(1, T);
time_it.CoD_R_GAP = cell(1, T);
time_it.CoD_Dome_GAP = cell(1, T);
time_it.CoD_Sta = cell(1, T);
time_it.CoD_G_Sta = cell(1, T);
time_it.CoD_Dome_Sta = cell(1, T);
time_it.CoD_Sta_GAP = cell(1, T);
time_it.CoD_G_Sta_GAP = cell(1, T);
% ɸѡʱ��
time_it_screen.CoD_GAP = cell(1,T);
time_it_screen.CoD_G_GAP = cell(1,T);
time_it_screen.CoD_R_GAP = cell(1,T);
time_it_screen.CoD_Dome_GAP = cell(1,T);
time_it_screen.CoD_Sta = cell(1,T);
time_it_screen.CoD_G_Sta = cell(1,T);
time_it_screen.CoD_Dome_Sta = cell(1,T);
time_it_screen.CoD_Sta_GAP = cell(1,T);
time_it_screen.CoD_G_Sta_GAP = cell(1,T);
% ɸѡ����
num_screen.CoD_GAP= cell(1, T);
num_screen.CoD_G_GAP= cell(1, T);
num_screen.CoD_R_GAP= cell(1, T);
num_screen.CoD_Dome_GAP= cell(1, T);
num_screen.CoD_Sta= cell(1, T);
num_screen.CoD_G_Sta= cell(1, T);
num_screen.CoD_Dome_Sta= cell(1, T);
num_screen.CoD_Sta_GAP= cell(1, T);
num_screen.CoD_G_Sta_GAP= cell(1, T);
%% =========================main loop======================================
for i = 1 : length(lambdas)
    lambda = lambdas(i);
    fprintf('\n ---- Regularization parameter %d / %d ----\n',i, length(lambdas))
    path_KL_solvers
end
save('result_MNIST-1-25', 'x_opt', 'x_it', 'primal_obj', 'theta_opt', 'stop_crit_it','Iter','time_it','time_it_screen','num_screen')














