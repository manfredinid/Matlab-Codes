function F=RR2008p_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,taurate,subsidyrate, cf, gcost)
% a_val is what HopenhaynRogerson1993 call n_{t-1}, aprime_val is n_t.

F=-Inf;

gnnlag=0; 
if aprime_val>a_val
    gnnlag=-gcost*((a_val-aprime_val)/a_val); % Note that gnnlag>=0, it is 'minus a negative number'
end
tau=taurate*(tau_val>=0)-subsidyrate*(tau_val<0);

% Physical capital:

% Labour
nbar=(((1-tau)*s_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));
    
pi =  p*(1-tau)*s_val*(aprime_val^alpha)*(nbar^gamma)-nbar-r*aprime_val-cf-gnnlag;

F=pi;

end