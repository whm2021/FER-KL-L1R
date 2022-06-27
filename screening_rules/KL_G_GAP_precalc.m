function  precalc= KL_G_GAP_precalc(A,y, lambda, precalc)
ind_y = find(y~=0);
Alpha = zeros(length(ind_y),1);
for i = 1 : length(ind_y)
    ind_i = ind_y(i);
    ind_a = A(ind_i,:)~=0;
    denominator_r = min(((lambda + precalc.NormA11(ind_a))./A(ind_i,ind_a)).^2);
    Alpha(i) = y(ind_i)./denominator_r;  
end
precalc.alpha = lambda^2 * min(Alpha);
