function F=capital_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma, delta, taurate,subsidyrate, cf, gcost)
% a_val is what HopenhaynRogerson1993 call n_{t-1}, aprime_val is n_t.

F=-Inf;

gnnlag=0;
if (1-delta)*aprime_val>a_val && a_val>0 % ROB: I ADDED THIS a_val>0
    gnnlag=(gcost/2)*(((1-delta)*aprime_val-a_val)/a_val)^2; % Note that gnnlag>=0, it is 'minus a negative number'
end
tau=taurate*(tau_val>0)-subsidyrate*(tau_val<0);

% Labour
nbar=((s_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));

pi =  p*s_val*(aprime_val^alpha)*(nbar^gamma)-nbar-r*(1-(-tau))*aprime_val-cf-gnnlag;

F=pi;

end