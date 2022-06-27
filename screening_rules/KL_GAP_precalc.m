function  precalc= KL_GAP_precalc(y, lambda, precalc)
denominator_r = (1 + max(precalc.NormA1, lambda).*precalc.NormAp1).^2;
denominator_r = denominator_r(y~=0);
precalc.alpha = lambda^2 * min(y(y~=0)./denominator_r);
