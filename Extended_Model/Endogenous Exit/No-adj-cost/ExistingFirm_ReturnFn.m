function F=ExistingFirm_ReturnFn(aprime_val, a_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma, delta,cf)

F=-Inf;


% Interest rate depends on earmarked vs non-earmarked credit firm
r=r_ear*psi_val + r_market*(1-psi_val);


kbar=(alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma*p)^(gamma/(1-gamma-alpha)) *(s_val)^(1/(1-alpha-gamma));


% Labour
nbar=(s_val*gamma*p)^(1/(1-gamma)) *kbar^(alpha/(1-gamma));


pi =  p*s_val*(kbar^alpha)*(nbar^gamma)-w*nbar-r*kbar-cf;

F=pi;

end