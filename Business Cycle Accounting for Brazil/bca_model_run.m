%This program computes a second-order  approximation to the policy functions of a simple neoclassical model (see ``Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function,'' by Stephanie Schmitt-Grohe and Martin Uribe, (JEDC, 2004, p. 755-775). The equilibrium conditions of the model can be written as: %E_t[f(yp,y,xp,x)=0, 

%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001, revised 22-Oct-2004

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = bca_model;

%Numerical Evaluation
%Steady State and Parameter Values

[b_y, delta, i_y, L, eta, RBAR, alph, nu, gama, n, A, TAUL, TAUK, TAUB, pssi, betta, phi, B,...
        c, cp, ik, ikp, la, lap, y, yp, k, kp, b, bp, r, rp, a, ap, tau_l, tau_lp, tau_k, tau_kp, tau_b, tau_bp]=bca_ss;


%Order of approximation desired 
approx = 2;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
