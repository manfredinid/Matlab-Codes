function [ys_, params, info] = model_observed_steadystate2(ys_, exo_, params)
% Steady state generated by Dynare preprocessor
    info = 0;
    ys_(1)=(params(7)*(exo_(2)*(1-params(1))/(params(1)*exo_(1)))^((1-params(1))*params(8))/(exo_(2)*(1/params(2)-1+params(6))/(params(1)*params(8))))^(1/(1-params(8)));
    ys_(4)=params(7)*ys_(1)^(params(1)*params(8))*(exo_(2)*(1-params(1))/(params(1)*exo_(1))*ys_(1))^((1-params(1))*params(8));
    ys_(3)=params(7)*ys_(1)^(params(1)*params(8))*(exo_(2)*(1-params(1))/(params(1)*exo_(1))*ys_(1))^((1-params(1))*params(8))-exo_(2)*(1-params(1))/(params(1)*exo_(1))*ys_(1)*exo_(1)*params(6)-ys_(1)*exo_(2)*params(6);
    ys_(2)=exo_(2)*(1-params(1))/(params(1)*exo_(1))*ys_(1);
    ys_(5)=ys_(1)+ys_(2);
    ys_(6)=ys_(4)-ys_(3);
    % Auxiliary equations
    check_=0;
end
