function [x, x_it, primal_obj_value, theta, stop_crit_it, iter, time_it, time_it_screen, num_screen] = CoD_KL_l1_Dome_GAP(A, y, a0, lambda, para, precalc)
% Author: Hongmei Wang
% Date: 24 Oct 2021

[n, m] = size(A);
%% ԭ����Ŀ�꺯��
%f.eval = @(a,b) sum(a.*log(a./b) - a + b); % KL distance
f.eval = @(a) sum(y(y~=0).*log(y(y~=0)./(a(y~=0)+ para.epsilon))) + sum(- y + a + para.epsilon); % force 0*log(0) = 0 (instead of NaN)
g.eval = @(a) lambda*norm(a, 1); % regularization
%% ��ʼ��
x0 = ones(m,1);
x = x0 + 0; % Signal to find (+0 avoid x to be just a pointer to x0)
ind_Fc=1:m; ind_F = [];

t = 1;
Ax = A*x;
primal_obj_value(1) =  f.eval(Ax) + g.eval(x) ;
stop_crit = inf;
stop_crit_it(1) = stop_crit;
x_it(:,1) = x0;
time_it(1) = 0;
time_it_screen(1)=0;

% CoD iterations
while (stop_crit > para.tol) && (t < para.max_iter)
    t = t + 1;
    fprintf('%4d,',t)
    x_old = x + 0; % +0 avoids x_old to be modified within the MEX function
    Ax = A(:,ind_Fc)*x(ind_Fc);
    tStart = tic;
    % ==========����ţ�ٷ�(Hsieh-Dhillon2011 Coord Descent)�������x========
    for i = 1:length(ind_Fc)
        t_coord = ind_Fc(i);
        tmp = y./(Ax+para.epsilon);
        tmp2 = tmp./(Ax+para.epsilon);
        g1 = A(:,t_coord).'*(1-tmp);
        g2 = (A(:,t_coord).^2).'*(tmp2);
        x(t_coord) = max(0, x(t_coord) -(g1+lambda)/g2);
        Ax = Ax + (x(t_coord)-x_old(t_coord))*A(:,t_coord);
    end
    %============================= ֹͣ׼�� ===============================
    % ����ԭ����Ŀ�꺯��ֵ
    primal = f.eval(Ax) + g.eval(x) ;
    % �����ż����Ŀ�꺯��ֵ
    theta = (y - Ax - para.epsilon)./(lambda*(Ax+para.epsilon));
    ATtheta = A.'*theta;
    theta = theta/max(1,max(ATtheta));
    dual = y(y~=0).'*log(1+lambda*theta(y~=0)) - sum(lambda*para.epsilon*theta);
    % �����ż��϶
    gap = primal - dual;
    % stop_crit = min(abs(gap),norm(x_old-x,2));
    stop_crit = gap;
    %============================ ɸѡ׼�� ================================
    tScreen = tic;
    center = theta;
    radius = sqrt(2*gap/precalc.alpha);
    s = (center'*a0 - 1)^2;
    num_s = 0;
    for j = 1 : length(ind_Fc)
        j_coord = ind_Fc(j);
        a = A(:,j_coord);
        if radius*a'*a0 > sqrt(s)*norm(a)
           pi = a'*center + radius*precalc.NormA(j_coord);
        else
            c = radius^2*norm(a0)^2 - s;
            AA = norm(a0)^2*c;
            BB = -2*(a'*a0)*c;
            CC = radius^2*(a'*a0)^2 - s*norm(a)^2;
            ss = (-BB + sqrt(BB^2-4*AA*CC))/(2*AA);
            
            aa0 = a - ss*a0;
            pi = aa0'*center + radius*norm(aa0(y~=0)) + ss;
        end
        if pi < 1-1e-8
            num_s = num_s + 1;
            x(j_coord) = 0;
            ind_F = [ind_F j_coord];
            
        end
    end
    time_it_screen(t) = toc(tScreen);
    
    ind_Fc = setdiff(ind_Fc, ind_F);
    num_screen(t) = num_s;
    %==================== �����Ҫ����Ľ�����и�ֵ =======================
    primal_obj_value(t) =  primal;
    x_it(:, t) = x;
    stop_crit_it(t) = stop_crit;
    iter = t-1;
    time_it(t) = toc(tStart); % total elapsed time since the beginning
end
end
