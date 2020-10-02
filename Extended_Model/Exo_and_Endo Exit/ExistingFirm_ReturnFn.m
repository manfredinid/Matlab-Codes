function F=ExistingFirm_ReturnFn(kprime_val, k_val,s_val, psi_val, tau_val, p,w,r_market,r_ear, alpha,gamma, delta, ctau, adjustcostparam)

F=-Inf;

adjcost=0;
if k_val>0 && kprime_val>k_val
    adjcost=(adjustcostparam/2)*(((kprime_val-k_val)/k_val)-delta)^2; % Note that adjcost>=0, it is 'minus a negative number'
end

% Interest rate depends on earmarked vs non-earmarked credit firm
r= tau_val*ctau + r_ear*psi_val + r_market*(1-psi_val);

% Labour
nbar=((s_val*p*gamma))^(1/(1-gamma)) *kprime_val^(alpha/(1-gamma));

pi =  p*s_val*(kprime_val^alpha)*(nbar^gamma)-w*nbar-r*kprime_val-adjcost-1;

F=pi;

end