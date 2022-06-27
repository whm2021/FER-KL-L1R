function [x, x_it, primal_obj_value, theta, stop_crit_it, iter, time_it, time_it_screen, num_screen] = CoD_KL_l1_Dome_Sta(A, y, a0, lambda, lambda0, x0, theta0, para, precalc)
% Author: Hongmei Wang
% Date: 17 Oct 2021
[n, m] = size(A);
%% 原问题目标函数
%f.eval = @(a,b) sum(a.*log(a./b) - a + b); % KL distance
f.eval = @(a) sum(y(y~=0).*log(y(y~=0)./(a(y~=0)+ para.epsilon))) + sum(- y + a + para.epsilon); % force 0*log(0) = 0 (instead of NaN)
g.eval = @(a) lambda*norm(a, 1); % regularization
%% 初始化
x = x0 + 0; % Signal to find (+0 avoid x to be just a pointer to x0)

t = 1;
%=============================== 筛选准则 =================================
tScreen = tic;
center = 0.5 * (1 + lambda/lambda0) * theta0;
g1 = y(y~=0).'*log(1+lambda*(lambda/lambda0)*theta0(y~=0)) - sum(lambda*para.epsilon*(lambda/lambda0)*theta0);
g2 = y(y~=0).'*log(1+lambda*theta0(y~=0)) - sum(lambda*para.epsilon*theta0);
g3 = lambda*(1-lambda/lambda0)*(para.epsilon-y./(1+lambda*theta0))'*theta0;
g4 = (1-lambda/lambda0)^2*norm(theta0,2)^2;
radius = sqrt(1/precalc.alpha*(-g1+g2+g3)-0.25*g4);

s = (center'*a0 - 1)^2;
for j = 1 : m
    a = A(:,j);
    if radius*a'*a0 > sqrt(s)*norm(a)
        pi(j) = a'*center + radius*precalc.NormA(j);
    else
        c = radius^2*norm(a0)^2 - s;
        AA = norm(a0)^2*c;
        BB = -2*(a'*a0)*c;
        CC = radius^2*(a'*a0)^2 - s*norm(a)^2;
        ss = (-BB + sqrt(BB^2-4*AA*CC))/(2*AA);
        
        aa0 = a - ss*a0;
        pi(j) = aa0'*center + radius*norm(aa0(y~=0)) + ss;
    end
end
ind_screen = (pi < 1);
x(ind_screen ~= 0) = 0;
time_it_screen(t) = toc(tScreen);
ind_Fc = find(ind_screen == 0);
num_screen = sum(ind_screen);
% ================================= end ===================================
Ax = A*x;
% 计算原问题目标函数值
primal = f.eval(Ax) + g.eval(x) ;
% 计算对偶问题目标函数值
theta = (y - Ax - para.epsilon)./(lambda*(Ax+para.epsilon));
ATtheta = A.'*theta;
theta = theta/max(1,max(ATtheta));
dual = y(y~=0).'*log(1+lambda*theta(y~=0)) - sum(lambda*para.epsilon*theta);
% 计算对偶间隙
gap = primal - dual;
%stop_crit = min(abs(gap),norm(x_old-x,2));
stop_crit = gap;

primal_obj_value(t) =  primal;
stop_crit_it(t) = stop_crit;
x_it(:,t) = x;
time_it(t) = 0;
iter = 0;
% CoD iterations
while (stop_crit > para.tol) && (t < para.max_iter)
    t = t + 1;
    fprintf('%4d,',t)
    x_old = x + 0; % +0 avoids x_old to be modified within the MEX function
    Ax = A(:,ind_Fc)*x(ind_Fc);
    tStart = tic;
    % ==========利用牛顿法(Hsieh-Dhillon2011 Coord Descent)逐个更新x========
    for i = 1:length(ind_Fc)
        t_coord = ind_Fc(i);
        tmp = y./(Ax+para.epsilon);
        tmp2 = tmp./(Ax+para.epsilon);
        g1 = A(:,t_coord).'*(1-tmp);
        g2 = (A(:,t_coord).^2).'*(tmp2);
        x(t_coord) = max(0, x(t_coord) -(g1+lambda)/g2);
        Ax = Ax + (x(t_coord)-x_old(t_coord))*A(:,t_coord);
    end
    %============================= 停止准则 ===============================
    % 计算原问题目标函数值
    primal = f.eval(Ax) + g.eval(x) ;
    % 计算对偶问题目标函数值
    theta = (y - Ax - para.epsilon)./(lambda*(Ax+para.epsilon));
    ATtheta = A.'*theta;
    theta = theta/max(1,max(ATtheta));
    dual = y(y~=0).'*log(1+lambda*theta(y~=0)) - sum(lambda*para.epsilon*theta);
    % 计算对偶间隙
    gap = primal - dual;
    %stop_crit = min(abs(gap),norm(x_old-x,2));
    stop_crit = gap;
    %==================== 对最后要保存的结果进行赋值 =======================
    primal_obj_value(t) =  primal;
    x_it(:, t) = x;
    stop_crit_it(t) = stop_crit;
    iter = t-1;
    time_it(t) = toc(tStart); % total elapsed time since the beginning
end
end