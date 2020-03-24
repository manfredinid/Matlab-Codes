function [tstar,sevec]=bca(x0,yt,ct,ht,it)

%  COPYRIGHT (c) 2003 BY PETER N. IRELAND.

% set starting values

rhoA = x0(1); % AR(1) Productivity
rhoL = x0(2); % AR(1) Labor Wedge
rhoK = x0(3); % AR(1) Capital Wedge
rhoB = x0(4); % AR(1) Bond Wedge
sigA = x0(5); % S.D. Productivity
sigL = x0(6); % S.D. Labor Wedge
sigK = x0(7); % S.D. Capital Wedge
sigB = x0(8); % S.D. Bond Wedge

% maximize likelihood

   bigtheto = [ rhoA , rhoL, rhoK, rhoB, sigA, sigL, sigK, sigB ];

  options = optimset('MaxIter',1000,'Display','Iter','MaxFunEvals',1000);

  thetstar = fminunc(@llfndbca,bigtheto,options,yt,ct,ht,it);
  
  % find standard errors

  thetstar = real(thetstar);

  rho1 = thetstar(1);
  rho2 = thetstar(2);
  rho3 = thetstar(3);
  rho4 = thetstar(4);
  sig1 = abs(thetstar(5));
  sig2 = abs(thetstar(6));
  sig3 = abs(thetstar(7));
  sig4 = abs(thetstar(8));

  tstar = [ rho1 , rho2, rho3, rho4, ...
             sig1, sig2, sig3, sig4]';

  scalinv=eye(8);

  tstars = tstar;
  
  fstar = llfndsebca(tstars,yt,ct,ht,it);

  eee = 1e-6;
  
  epsmat = eee*eye(8);

  hessvec = zeros(8,1);

  for i = 1:8

    hessvec(i) = llfndsebca(tstars+epsmat(:,i),yt,ct,ht,it);

  end

  hessmat = zeros(8,8);

  for i = 1:8

    for j = 1:8

      hessmat(i,j) = (llfndsebca(tstars+epsmat(:,i)+epsmat(:,j),yt,ct,ht,it) ...
                        -hessvec(i)-hessvec(j)+fstar)/eee^2;

    end

  end

  bighx = scalinv*inv(hessmat)*scalinv';

  sevec = sqrt(diag(bighx));
