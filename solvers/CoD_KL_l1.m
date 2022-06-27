function [x, x_it, primal_obj_value, theta, stop_crit_it, iter, time_it] = CoD_KL_l1(A, y, lambda, para)
% Author: Hongmei Wang
% Date: 13 Oct 2021
[n, m] = size(A);
%% 原问题目标函数
%f.eval = @(a,b) sum(a.*log(a./b) - a + b); % KL distance
f.eval = @(a) sum(y(y~=0).*log(y(y~=0)./(a(y~=0)+ para.epsilon))) + sum(- y + a + para.epsilon); % force 0*log(0) = 0 (instead of NaN)
g.eval = @(a) lambda*norm(a, 1); % regularization
%% 初始化
x0 = ones(m,1);
x = x0 + 0; % Signal to find (+0 avoid x to be just a pointer to x0)

t = 1;
Ax = A*x;
primal_obj_value(1) =  f.eval(Ax) + g.eval(x) ;
stop_crit = inf;
stop_crit_it(1) = stop_crit;
x_it(:,1) = x0;
time_it(1) = 0;

% CoD iterations
while (stop_crit > para.tol) && (t < para.max_iter)
    t = t + 1;
    fprintf('%4d,',t) 
    x_old = x + 0; % +0 avoids x_old to be modified within the MEX function
    tStart = tic;
    % ==========利用牛顿法(Hsieh-Dhillon2011 Coord Descent)逐个更新x======== 
    for t_coord = 1:m
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
