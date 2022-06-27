function  alpha= KL_R_GAP_precalc(y, lambda, theta, radiu)
ind_y = find(y~=0);
denominator_r =(1 + lambda.*(theta(ind_y) + radiu)).^2;
alpha = lambda^2 *min(y(ind_y)./denominator_r);