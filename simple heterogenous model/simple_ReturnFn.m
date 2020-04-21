function F=simple_ReturnFn(aprime_val, a_val, s_val, p, alpha, cf)
% Note that neither aprime_val nor a_val is actually used for anything in
% this model. But VFI toolkit is not set up to handle that they do not exist.

F=-Inf;

kbar = (s_val*p*alpha)^(1/(1-alpha));

% This example is taken from Chris Edmonds lecture notes.
pi=p*s_val*(kbar^alpha)-kbar-cf;

F=pi;

end