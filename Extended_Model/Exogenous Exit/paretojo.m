function [ergoeps, pie] = paretojo(neps, eps, shape, rhoeps)

% This computes the ergodic distribution and the transition matrix of a
% bounded Pareto distribution, given grids of epsilon values

ergoeps = zeros(1,neps);
pie = zeros(neps,neps);

% emin and emax are imaginary end points of a continuous distribution
emin = eps(1) - (eps(2)-eps(1))/2;
emax = eps(neps) + (eps(neps)-eps(neps-1))/2;

mideps = zeros(1,neps-1);
for i = 1:neps-1
    mideps(i) = (eps(i)+eps(i+1))/2;
end

% compute density or mass at each epsilon point
ergoeps(1) = bddparetocdf(emin, emax, shape, mideps(1));
for i = 2:neps-1
    mass = bddparetocdf(emin,emax,shape, mideps(i-1));
    ergoeps(i) = bddparetocdf(emin,emax,shape, mideps(i)) - mass;
end
ergoeps(neps) = 1.0 - bddparetocdf(emin, emax, shape, mideps(neps-1));


% transition matrix given rhoeps
for i = 1:neps
    pie(i,i) = rhoeps;
    pie(i,:) = pie(i,:) + (1.0-rhoeps)*ergoeps;
end








