function [at,tht,tkt,tbt]=bca_wedge_alternative2(x1,gx,hx,yt,ct,ht,it,kt)

% This program computes the wedges of the model

% Parameters of the model

[b_y, delta, i_y, L, eta, RBAR, alph, nu, gama, n, A, TAUL, TAUK, TAUB, pssi, betta, phi, B,...
        c, cp, ik, ikp, la, lap, y, yp, k, kp, b, bp, r, rp, a, ap, tau_l, tau_lp, tau_k, tau_kp, tau_b, tau_bp]=bca_ss;

% Initial Series 

T=length(yt);

at=zeros(T,1);

tht=zeros(T,1);

tkt=zeros(T,1);

tbt=zeros(T,1);

xt = [zeros(3,T); at'; tht'; tkt'; tbt'];

xt = [xt, zeros(7,1)];

ythat = log(yt) - log(mean(yt));

ithat = log(it) - log(mean(it));

hthat = log(ht) - log(mean(ht));

cthat = log(ct) - log(mean(ct));
 
zt = [cthat, ithat, hthat, ythat]';

for kk=1:T
    
    B = [zt(:,kk)] - [gx(:,1:3) * xt(1:3,kk)];
    
    A = [gx(:,4:7)];
    
    vx = A\B;

    at(kk) = vx(1);
    
    tht(kk) = vx(2);
    
    tkt(kk) = vx(3);
    
    tbt(kk) = vx(4);
    
    xt(4:7,kk) = [at(kk); tht(kk); tkt(kk); tbt(kk)];
    
    xtp = hx * xt(:,kk);
    
    xt(1:3,kk+1) = xtp(1:3);
    
end
