function [residual, g1, g2, g3] = model_observed_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 6, 1);

%
% Model equations
%

psi__ = x(2)*(1-params(1))/(params(1)*x(1));
A__ = params(4)^params(1)*params(3)^(1-params(1));
expH__ = params(1)*params(8);
expL__ = (1-params(1))*params(8);
T28 = A__*y(1)^expH__;
T32 = (y(1)*psi__)^expL__;
T33 = T28*T32;
T38 = y(3)^(-params(5));
T51 = expH__*T33/(y(1)*(x(2)+x(1)))+expL__*T33/(y(1)*psi__*(x(2)+x(1)))+1-params(6);
lhs =y(4);
rhs =T33;
residual(1)= lhs-rhs;
lhs =T38;
rhs =T38*params(2)*T51;
residual(2)= lhs-rhs;
lhs =y(3);
rhs =T33-x(1)*(y(1)*psi__-y(1)*psi__*(1-params(6)))-x(2)*(y(1)-y(1)*(1-params(6)));
residual(3)= lhs-rhs;
lhs =y(2);
rhs =x(2)*(1-params(1))/(params(1)*x(1))*y(1);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(1)+y(2);
residual(5)= lhs-rhs;
lhs =y(6);
rhs =y(4)-y(3);
residual(6)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(6, 6);

  %
  % Jacobian matrix
  %

T78 = T32*A__*getPowerDeriv(y(1),expH__,1)+T28*psi__*getPowerDeriv(y(1)*psi__,expL__,1);
T105 = getPowerDeriv(y(3),(-params(5)),1);
  g1(1,1)=(-T78);
  g1(1,4)=1;
  g1(2,1)=(-(T38*params(2)*((y(1)*(x(2)+x(1))*expH__*T78-expH__*T33*(x(2)+x(1)))/(y(1)*(x(2)+x(1))*y(1)*(x(2)+x(1)))+(y(1)*psi__*(x(2)+x(1))*expL__*T78-expL__*T33*psi__*(x(2)+x(1)))/(y(1)*psi__*(x(2)+x(1))*y(1)*psi__*(x(2)+x(1))))));
  g1(2,3)=T105-T51*params(2)*T105;
  g1(3,1)=(-(T78-x(1)*(psi__-psi__*(1-params(6)))-x(2)*(1-(1-params(6)))));
  g1(3,3)=1;
  g1(4,1)=(-(x(2)*(1-params(1))/(params(1)*x(1))));
  g1(4,2)=1;
  g1(5,1)=(-1);
  g1(5,2)=(-1);
  g1(5,5)=1;
  g1(6,3)=1;
  g1(6,4)=(-1);
  g1(6,6)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,36);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,216);
end
end
end
end
