function llfn = llfndbca(bigthet,yt,ct,ht,it);

%  COPYRIGHT (c) 2003 BY PETER N. IRELAND.

% define variables and parameters

bigthet = real(bigthet);

rhoA = bigthet(1);
rhoL = bigthet(2);
rhoK = bigthet(3);
rhoB = bigthet(4);
sigA = bigthet(5);
sigL = bigthet(6);
sigK = bigthet(7);
sigB = bigthet(8);
 
% untransform parameters

rhoA = rhoA;
rhoL = rhoL;
rhoK = rhoK;
rhoB = rhoB;
sigA = abs(sigA);
sigL = abs(sigL);
sigK = abs(sigK);
sigB = abs(sigB);


% Solve model

bca_model_run;

% form matrices FX, GX, and QX

  bigfx = hx;

  biggx = gx;
  
  bigbx = [zeros(3,4); eye(4,4)];

  bigv1x = [sigA^2, 0, 0, 0;
               0,sigL^2,0,0;
               0,0,sigK^2,0;
               0,0,0,sigB^2];

  bigqx = [bigbx*bigv1x*bigbx'];

% put data in deviation form

  bigt = length(yt);

  ythat = log(yt) - log(mean(yt));

  ithat = log(it) - log(mean(it));

  hthat = log(ht) - log(mean(ht));

  cthat = log(ct) - log(mean(ct));
 
  dthat = [cthat, ithat, hthat, ythat];

% evaluate negative log likelihood

  xt = zeros(7,1);

  bigsig1 = inv(eye(49)-kron(bigfx,bigfx))*bigqx(:);

  bigsigt = reshape(bigsig1,7,7);

  llfn = (3*bigt/2)*log(2*pi);

  for t = 1:bigt

    ut = dthat(t,:)' - biggx*xt;

    omegt = biggx*bigsigt*biggx';

    omeginvt = inv(omegt);

    llfn = llfn + (1/2)*(log(det(omegt))+ut'*omeginvt*ut);

    bigkt = bigfx*bigsigt*biggx'*omeginvt;

    xt = bigfx*xt + bigkt*ut;

    bigsigt = bigqx + bigfx*bigsigt*bigfx' ...
                - bigfx*bigsigt*biggx'*omeginvt*biggx*bigsigt*bigfx';

  end

% penalize constraint violations

  if abs(imag(llfn)) > 0

    llfn = real(llfn) + 1e8;

  end

  if max(abs(eig(hx))) >= 1

    llfn = llfn + 1e8;

  end
  