function [residual, g1, g2, g3] = model_observed_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(6, 1);
psi__ = x(it_, 2)*(1-params(1))/(params(1)*x(it_, 1));
A__ = params(4)^params(1)*params(3)^(1-params(1));
expH__ = params(1)*params(8);
expL__ = (1-params(1))*params(8);
T28 = A__*y(1)^expH__;
T32 = (y(1)*psi__)^expL__;
T33 = T28*T32;
T42 = params(2)*y(8)^(-params(5));
T45 = A__*y(2)^expH__;
T47 = (psi__*y(2))^expL__;
T48 = T45*T47;
T59 = expH__*T48/(y(2)*(x(it_, 2)+x(it_, 1)))+expL__*T48/(psi__*y(2)*(x(it_, 2)+x(it_, 1)))+1-params(6);
lhs =y(5);
rhs =T33;
residual(1)= lhs-rhs;
lhs =y(4)^(-params(5));
rhs =T42*T59;
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T33-x(it_, 1)*(psi__*y(2)-y(1)*psi__*(1-params(6)))-x(it_, 2)*(y(2)-y(1)*(1-params(6)));
residual(3)= lhs-rhs;
lhs =y(3);
rhs =x(it_, 2)*(1-params(1))/(params(1)*x(it_, 1))*y(2);
residual(4)= lhs-rhs;
lhs =y(6);
rhs =y(2)+y(3);
residual(5)= lhs-rhs;
lhs =y(7);
rhs =y(5)-y(4);
residual(6)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(6, 10);

  %
  % Jacobian matrix
  %

T82 = getPowerDeriv(y(1)*psi__,expL__,1);
T86 = T32*A__*getPowerDeriv(y(1),expH__,1)+T28*psi__*T82;
T98 = getPowerDeriv(psi__*y(2),expL__,1);
T102 = T47*A__*getPowerDeriv(y(2),expH__,1)+T45*psi__*T98;
T135 = T28*T82*y(1)*(-(params(1)*x(it_, 2)*(1-params(1))))/(params(1)*x(it_, 1)*params(1)*x(it_, 1));
T137 = y(2)*(-(params(1)*x(it_, 2)*(1-params(1))))/(params(1)*x(it_, 1)*params(1)*x(it_, 1));
T139 = T45*T98*T137;
T167 = y(2)*(1-params(1))/(params(1)*x(it_, 1));
  g1(1,1)=(-T86);
  g1(1,5)=1;
  g1(1,9)=(-T135);
  g1(1,10)=(-(T28*T82*y(1)*(1-params(1))/(params(1)*x(it_, 1))));
  g1(2,2)=(-(T42*((y(2)*(x(it_, 2)+x(it_, 1))*expH__*T102-expH__*T48*(x(it_, 2)+x(it_, 1)))/(y(2)*(x(it_, 2)+x(it_, 1))*y(2)*(x(it_, 2)+x(it_, 1)))+(psi__*y(2)*(x(it_, 2)+x(it_, 1))*expL__*T102-expL__*T48*psi__*(x(it_, 2)+x(it_, 1)))/(psi__*y(2)*(x(it_, 2)+x(it_, 1))*psi__*y(2)*(x(it_, 2)+x(it_, 1))))));
  g1(2,4)=getPowerDeriv(y(4),(-params(5)),1);
  g1(2,8)=(-(T59*params(2)*getPowerDeriv(y(8),(-params(5)),1)));
  g1(2,9)=(-(T42*((y(2)*(x(it_, 2)+x(it_, 1))*expH__*T139-y(2)*expH__*T48)/(y(2)*(x(it_, 2)+x(it_, 1))*y(2)*(x(it_, 2)+x(it_, 1)))+(psi__*y(2)*(x(it_, 2)+x(it_, 1))*expL__*T139-expL__*T48*(psi__*y(2)+(x(it_, 2)+x(it_, 1))*T137))/(psi__*y(2)*(x(it_, 2)+x(it_, 1))*psi__*y(2)*(x(it_, 2)+x(it_, 1))))));
  g1(2,10)=(-(T42*((y(2)*(x(it_, 2)+x(it_, 1))*expH__*T45*T98*T167-y(2)*expH__*T48)/(y(2)*(x(it_, 2)+x(it_, 1))*y(2)*(x(it_, 2)+x(it_, 1)))+(psi__*y(2)*(x(it_, 2)+x(it_, 1))*expL__*T45*T98*T167-expL__*T48*(psi__*y(2)+(x(it_, 2)+x(it_, 1))*T167))/(psi__*y(2)*(x(it_, 2)+x(it_, 1))*psi__*y(2)*(x(it_, 2)+x(it_, 1))))));
  g1(3,1)=(-(T86-x(it_, 1)*(-(psi__*(1-params(6))))-x(it_, 2)*(-(1-params(6)))));
  g1(3,2)=(-((-(x(it_, 1)*psi__))-x(it_, 2)));
  g1(3,4)=1;
  g1(3,9)=(-(T135-(psi__*y(2)-y(1)*psi__*(1-params(6))+x(it_, 1)*(T137-(1-params(6))*y(1)*(-(params(1)*x(it_, 2)*(1-params(1))))/(params(1)*x(it_, 1)*params(1)*x(it_, 1))))));
  g1(3,10)=(-(T28*T82*y(1)*(1-params(1))/(params(1)*x(it_, 1))-x(it_, 1)*(T167-(1-params(6))*y(1)*(1-params(1))/(params(1)*x(it_, 1)))-(y(2)-y(1)*(1-params(6)))));
  g1(4,2)=(-(x(it_, 2)*(1-params(1))/(params(1)*x(it_, 1))));
  g1(4,3)=1;
  g1(4,9)=(-T137);
  g1(4,10)=(-T167);
  g1(5,2)=(-1);
  g1(5,3)=(-1);
  g1(5,6)=1;
  g1(6,4)=1;
  g1(6,5)=(-1);
  g1(6,7)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],6,100);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],6,1000);
end
end
end
end
