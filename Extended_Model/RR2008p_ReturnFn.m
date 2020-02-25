function F=HopenhaynRogerson1993_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,taurate,subsidyrate, cf)
% a_val is what HopenhaynRogerson1993 call n_{t-1}, aprime_val is n_t.

F=-Inf;

tau=taurate*(tau_val>=0)-subsidyrate*(tau_val<0);

% Physical capital:
kbar=(alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma)^(gamma/(1-gamma-alpha)) *(p*s_val*(1-tau))^(1/(1-alpha-gamma));

% Labour
nbar=(((1-tau)*s_val*p*gamma))^(1/(1-gamma)) *kbar^(alpha/(1-gamma));
    
pi =  p*(1-tau)*s_val*(kbar^alpha)*(nbar^gamma)-nbar-r*kbar-p*cf;

F=pi;

end