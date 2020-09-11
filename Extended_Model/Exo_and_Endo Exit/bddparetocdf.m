function cdf = bddparetocdf(emin, emax, shape, eval)

% This calculates the cumulative density (or mass) at eval for a
% discretized bounded Pareto distribution over (emin,emax)

paretocdf = 1.0 - (emin/eval)^shape;
cdf = paretocdf/(1.0 - (emin/emax)^shape);

end