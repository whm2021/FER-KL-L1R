clear;clear all
load('result_test.mat')
load('test_data.mat')
%% 计算总时间
[~, n] = size(A);
m=100;
t = size(1,m);
% CoD
for i = 1:m
    t(i) = sum(time_it.CoD{1,i});   
end
final_time_all.CoD = sum(t);
% CoD_GAP
for i = 1:m
    t(i) = sum(time_it.CoD_GAP{1,i});   
end
final_time_all.CoD_GAP = sum(t);
% CoD_G_GAP
for i = 1:m
    t(i) = sum(time_it.CoD_G_GAP{1,i});   
end
final_time_all.CoD_G_GAP = sum(t);
% CoD_R_GAP
for i = 1:m
    t(i) = sum(time_it.CoD_R_GAP{1,i});   
end
final_time_all.CoD_R_GAP = sum(t);
% CoD_G_Sta
for i = 1:m
    t(i) = sum(time_it.CoD_G_Sta{1,i});   
end
final_time_all.CoD_G_Sta = sum(t);
% CoD_G_Sta_GAP
for i = 1:m
    t(i) = sum(time_it.CoD_G_Sta_GAP{1,i});   
end
final_time_all.CoD_G_Sta_GAP = sum(t);
%% Speedup
% CoD
% CoD_GAP
final_speedup.CoD_GAP = final_time_all.CoD/final_time_all.CoD_GAP;
% CoD_G_GAP
final_speedup.CoD_G_GAP = final_time_all.CoD/final_time_all.CoD_G_GAP;
% CoD_R_GAP
final_speedup.CoD_R_GAP = final_time_all.CoD/final_time_all.CoD_R_GAP;
% CoD_G_Sta
final_speedup.CoD_G_STA = final_time_all.CoD/final_time_all.CoD_G_Sta;
% CoD_G_Sta_GAP
final_speedup.CoD_G_Sta_GAP = final_time_all.CoD/final_time_all.CoD_G_Sta_GAP;
%% 筛选比例
m=100;
n_s = size(1,m);
% CoD_GAP
for i = 1:m
    n_s_GAP(i) = sum(num_screen.CoD_GAP{1,i});   
end
final_num_screen.CoD_GAP = sum(n_s_GAP)/(n*100);
% CoD_G_GAP
for i = 1:m
    n_s_G_GAP(i) = sum(num_screen.CoD_G_GAP{1,i});   
end
final_num_screen.CoD_G_GAP = sum(n_s_G_GAP)/(n*100);
% CoD_R_GAP
for i = 1:m
    n_s_R_GAP(i) = sum(num_screen.CoD_R_GAP{1,i});   
end
final_num_screen.CoD_R_GAP = sum(n_s_R_GAP)/(n*100);
% CoD_G_Sta
n_s_G_STA(1) = 0;
for i = 2:m
    n_s_G_STA(i) = num_screen.CoD_G_Sta{1,i};   
end
final_num_screen.CoD_G_STA = sum(n_s_G_STA)/(n*100);
% CoD_G_Sta_GAP
n_s_G_Sta_GAP(1) = sum(num_screen.CoD_G_Sta_GAP{1,1});
for i = 2:m
    n_s_G_Sta_GAP_sta(i) = num_screen.CoD_G_Sta_GAP{1,i}.Sta;
    n_s_G_Sta_GAP_gap(i) = sum(num_screen.CoD_G_Sta_GAP{1,i}.GAP); 
    n_s_G_Sta_GAP(i) = n_s_G_Sta_GAP_sta(i) + n_s_G_Sta_GAP_gap(i);
end
final_num_screen.CoD_G_Sta_GAP = sum(n_s_G_Sta_GAP)/(n*100);
%% 计算原问题最优解之间的误差
m=100;
% CoD_GAP
for i = 1 : m
    KL(:,i) = x_opt.CoD{1,i};
    GAP(:,i) = x_opt.CoD_GAP{1,i};
end
final_error.CoD_GAP = norm(KL-GAP);
% CoD_G_GAP
for i = 1 : m
    G_GAP(:,i) = x_opt.CoD_G_GAP{1,i};
end
final_error.CoD_G_GAP = norm(KL-G_GAP);
% CoD_R_GAP
for i = 1 : m
    R_GAP(:,i) = x_opt.CoD_R_GAP{1,i};
end
final_error.CoD_R_GAP = norm(KL-R_GAP);
% CoD_G_Sta
for i = 1 : m
    G_STA(:,i) = x_opt.CoD_G_Sta{1,i};
end
final_error.CoD_G_STA = norm(KL-G_STA);
% CoD_G_Sta_GAP
for i = 1 : m
    FER(:,i) = x_opt.CoD_G_Sta_GAP{1,i};
end
final_error.FER = norm(KL-FER);
