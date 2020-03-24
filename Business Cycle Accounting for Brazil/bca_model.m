function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = bca_model

%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001, revised 22-Oct-2004

%Define parameters
syms b_y delta i_y L eta RBAR alph nu rhoA rhoL rhoK rhoB gama n A TAUL TAUK TAUB pssi betta phi B

%Define variables 
syms c cp ik ikp la lap y yp k kp b bp r rp a ap tau_l tau_lp tau_k tau_kp tau_b tau_bp 

%Give functional form for  production, and utility functions

%Utility function 
uf = log (c) + pssi * log (1 - la);
ufp = log (cp) + pssi * log (1 - lap);

%Production Function 
pf = a * ( k ^ alph ) * ( la ^ ( 1 - alph ) );
pfp = ap * ( kp ^ alph ) * ( lap ^ ( 1 - alph ) );

%Adjustment Costs
phif = phi / 2 * ( ik / k - delta - gama - n - gama * n ) ^ 2;
phifp = phi / 2 * ( ikp / kp - delta - gama - n - gama * n ) ^ 2;

%Partial Derivatives

%Marginal Utility 
muc = diff(uf,'c'); 
mucp = diff(ufp,'cp'); 

%Marginal Disutility Labor
mula = diff(uf,'la'); 
mulap = diff(ufp,'lap'); 

%Marginal Productivity Labor
mpla = diff(pf,'la'); 
mplap = diff(pfp,'lap'); 

%Marginal Productivity of Capital
mpk = diff(pf,'k'); 
mpkp = diff(pfp,'kp'); 

%Marginal Adjustment Cost
mphif = phi * ( ik / k - delta - gama - n - gama * n );
mphifp = phi * ( ikp / kp - delta - gama - n - gama * n );

%Write equations fi, i=1:11

f1 = mula / muc + mpla * tau_l; % Consumption Leisure Decision

f2 = muc / ( 1 - mphif ) - ( betta / (1+gama) ) * ( mucp * ( tau_kp * mpkp + ...
    1/ ( 1 - mphifp ) * ( ( 1 - delta ) - phifp + mphifp * ikp / kp ) ) ); % Capital Euler Equation
f3 = muc - ( betta / (1 + gama) ) * ( mucp * tau_bp * rp ); % Euler Equation for bond
f4 = (1 + gama) * (1 + n) * bp - ( c + ik + b * r - pf ); % Accumulation of Net Foreign Assets
f5 = (1 + gama) * (1 + n) * kp - ( (1 - delta) * k + ik - phif * k ); % Law of Motion for Capital


f6 = rp - RBAR * (bp/B) ^ eta; % Supply of funds

f7 = y - pf; % Definition of Output


f8 = log(ap) - rhoA * log(a); % Productivity
f9 = log(tau_lp) - rhoL * log(tau_l); % Labor Wedge
f10 = log(tau_kp) - rhoK * log(tau_k); % Capital Wedge
f11 = log(tau_bp/TAUB) - rhoB * log(tau_b/TAUB); % Bond Wedge

%Create function f

f = [f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11];

% c cp ik ikp la lap y yp k kp b bp r rp a ap tau_l tau_lp tau_k tau_kp tau_b tau_bp 

% Define the vector of controls, y, and states, x
x = [b k r a tau_l tau_k tau_b];
y = [c ik la y];
xp = [bp kp rp ap tau_lp tau_kp tau_bp];
yp = [cp ikp lap yp];

%Make f a function of the logarithm of the state and control vector
f = subs(f, [x,y,xp,yp], exp([x,y,xp,yp]));

%Compute analytical derivatives of f
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);