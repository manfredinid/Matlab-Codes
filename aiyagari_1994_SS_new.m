%%%---------------
% Aiyagari(1994) - replication code
% Heterogeneous agent model (no aggregate shocks)
% heterogeneous households


% from Hugget to Aiygari: add labor income process and aggregate production
% function

% production model where agents differ in that they face uninsured 
% idiosyncratic labor endowment shocks and trade assets between themselves

% Code reference 

%% Model
%Households (of measure 1):
%1. Own capital
%2. Supply labor (exogenous and stochastic)
%3. Consume
% Firms: Rent capital and hire labor to produce
% Prices (r and w) are taken as given by households and firms

% ex-ante the households are identical, but they are hit by different
% shocks that change the individual decisions
%%%%--------------
clear all;
close all;
tic;
%% calibrated model parameters
%
sigma  = 2;            % relative risk aversion (curvature of utility function)      
beta   = 0.9;          % subjective discount factor 
delta  = 0.04;         % depreciation rate
alpha  = 0.3;         % capital's share of income



%% Grid features

%  productivity shock - labor states
% states
theta  = [1 0.01];            % Exogenous states
nlap = length(theta);

%Transition probabilities
% for the idiosyncratic productivity shock that follows a stochastic Markov process
% Markov chain
prob   = [ 2/3 1/3; 1/3 2/3]; % prob(i,j) = probability (s(t+1)=sj | s(t) = si)

%check matrix
[ev,ed] = eig(prob);
[emax,inmax] = max(diag(ed));
if emax~=1
   disp('are you sure the matrix prob is correct?');
end

%   form capital grid
%   
maxkap = 5;                      % maximum value of capital grid   
inckap = 5/50;                   % size of capital grid increments
nkap   = 50;                     % number of grid points
kap=linspace(0,maxkap,nkap);     % state of assets 

%% Finding a Steady State Equilibrium
% Matlab. We concentrate on stationary equilibria in the sense that we 
% require the interest rate (and therefore also the labor wage, wt = w) to be
% constant rt = r and consistent with capital market clearing.

% 1. set the tolerance parameter epsilon
% 2. guess the aggregate interest rate
% 3. compute the capital demand and the wage using the FOC
% 4. solve the problem of the agent and obtaion the optitmal decision rule
% for capital stock holding

% 5. using the political function and the law of montion for s find the
% distribution x(s,k)
% 6. compute the interest rate
% 7. Compare r0 and r1 until |r0 - r1| < epsilon



%% STEPS 1 TO 4
% In step 4
% Derive the household’s optimality conditions with respect to consumption,
% labor supply, and future asset holdings, taking as given factor prices, 
% wt and rt, and the transition probabilities.

% to solve step 4
% 1. set an arbitraty upperbond for the space of capital to make the domain of
% the VFI compact
% 2. approximate the value fucntion 
% 3. set the tolerance parameter epsilon
% 4. set an inicial guess for the value function
% 5. solve the Bellman equation
% 6. compare the inicial guess and the update value


liter   = 0;
maxiter = 50;
toler   = 1e-6;
metric  = 10;


r_L=-delta;
r_H = beta^(-1) - 1;


% Labor

D = zeros(length(prob));
[ev,ed] = eig(prob);
[emax,inmax] = max(diag(ed));
D(inmax,inmax) = emax;
pinf = ev*D*inv(ev);
L = theta*max(pinf)';

%%
while  (metric > toler) && (liter <= maxiter)
  
% Guess on r* (equilibrium interest rate) 
r = (r_L+r_H)/2;    


%  calculate rental rate of capital and wage
% from the firm's optimization problem
% k is K/L
k = ((r+delta)/alpha)^(1/(alpha-1));
w = (1-alpha)*k^(alpha);
util=zeros(nlap,nkap,nkap); 
   

   Fobj=zeros(nlap, nkap, nkap);
   v = zeros(nlap, nkap);
   Tv = zeros(nlap, nkap);
   Tg = zeros(nlap, nkap);
   dif=10;
   
 
   for ii=1:nlap
  for i=1:nkap
         kap1=(i-1)*inckap;
         for j=1:nkap
             kapp = (j-1)*inckap;
               c = theta(ii)*w + (r + (1-delta))*kap1 - kapp; 
               if c > 0
                  util(ii,j,i)=(c)^(1-sigma)/(1-sigma);
               else
            %  tabulate the utility function such that for zero or negative
            %  consumption utility remains a large negative number so that
            %  such values will never be chosen as utility maximizing 
                   util(ii,j,i)= - Inf;
               end
         end
  end
   end


      
   while dif ~= 0;   
   for ii=1:nlap
    for i=1:nkap
        for j=1:nkap
           Fobj(ii,j,i)=util(ii,j,i)+beta*(prob(ii,:)*v(:,i));
       its = 0+i;
        end
    end
   end
   
   for ii = 1:nlap
       for i = 1:nkap
           % maximum value and value of x at which this maximum is attained
    [Tv(ii,i), Tg(ii,i)]= max(Fobj(ii, i, :));
       end
   end

       dif =max(any(Tv- v));
       v=Tv;
       decis = Tg;
       
       
   end  
   decis = (decis-1)*inckap;
  
       
 
   %% STEPS 1 TO 7
   % 5. using the political function and the law of montion for s find the
   % distribution x(s,k)
   % 6. compute the interest rate
   % 7. Compare r0 and r1 until |r0 - r1| < epsilon
   
   
 
   Q = zeros(nlap*nkap);
  for i = 1:nlap
      for ii = 1:nkap
          for j = 1:nlap
              for jj = 1:nkap
       Q((i-1)*nkap+ii,(j-1)*nkap+jj)=0;
       if jj == Tg(i,ii)
                      Q((i-1)*nkap+ii,(j-1)*nkap+jj) = prob(i,j);
       end
              end
          end 
      end
  end
      
  
  
  n = length(Q);
  Pi1=(1/(n))*ones(1,n);
  dif2=10;
  while dif2>toler
      pi_ss = Pi1*Q;
      dif2=(max(abs(pi_ss-Pi1)));
      Pi1=pi_ss;
  end
  Pi_ss = Pi1';
  
 kk=decis(:);
 meanK= Pi_ss'*kk;
 
   %  calculate measure over (k,s) pairs
   %  lambda has same dimensions as decis
   %
   lambda=zeros(nkap,2);
   lambda(:)=Pi_ss;
   
    [v1,d1]=eig(prob');
   [dmax,imax]=max(diag(d1));
   probst1=v1(:,imax);
   ss=sum(pi_ss);
   probst1=pi_ss/ss;
   probk=sum(lambda');     %  stationary distribution of `captal' 
   probk=probk'; 
   
   
    
    %%

K = k*L;   

d = K - meanK;

if d==0

elseif d>0
		r_L = r;
elseif d<0
		r_H = r;
end


metric = abs(r_H - r_L);

liter = liter+1;

disp('iteration')

disp(liter)

K= full(K);
complete_market =((( beta^(-1) - 1)+delta)/alpha)^(1/(alpha-1)*L);
end

%%
%   print out results
%
disp('PARAMETER VALUES');
disp('');
disp('    sigma      beta      delta      alpha      '); 
disp([     sigma       beta     delta           alpha ]);
disp(''); 
disp('EQUILIBRIUM RESULTS ');
disp('');
disp('  interest  capital-labor   wage ');
disp([     r                       k                       w]);
disp('  capital      labor');
disp([     K           L]);
disp('  precautionary capital stock  ');
disp([ K-complete_market])

disp('the end')
figure;
bar(1:length(Pi_ss), Pi_ss, 'k')
title('Distribution of individuals')

 Pi_0 = Pi_ss(1:50);
 Pi_1 = Pi_ss(51:100);
 figure;
 bar(1:length(Pi_0),  Pi_0+Pi_1, 'k')
 title('Distribution of Capital')
 ylabel('% of agents');
 
 figure;
plot(kap,probk);
title('Distribution of Capital');
ylabel('% of agents');

 
 
 g0 = Tg(1,:); 
 g1 = Tg(2,:);
 
 figure;
plot(kap, kap(g0), 'k')
hold on
plot(kap, kap(g1), 'blue')
hold on
plot(kap, kap, '-r' )
legend('theta = 0.01', 'theta=1', 'reference','Location','SouthEast')
title('Policy Function')

toc;


%% One Agent Simulation
disp('ONE AGENT SIMULATION')
%Start from any initial condition for one individual
% Simulate they behavior of this individual for T periods, say 105,000
%Aggregate over this one guy’s policy function from t = 5001 to T. 
%This should approximate the stationary distribution
% If you want to check whether this works, simulate the same guy for T more periods. 
%Check if the to distributions are the same
clear iter capital_sim state_sim sz sZ rp N T
N=20;
T =20;  % number of periods to simulate

% randon number from a uniform distribution
rand('state',1); % reset the random number generator
rand_process=rand(N,T);

% Initialization matrix
state_sim    = zeros(N,T) ;
capital_sim    = zeros(N,T) ;
sZ      = ones(N,T) ;


if N>1
    sz(:,1) = [ 1*ones(round(N*1/2),1)
                2*ones(round(N*1/2),1) ];
else
     sZ(:,1) = 1 ;
end

state_sim(:,1) = theta(sZ(1)) ;

capital_sim(:,1) = K-inckap;


rand('state',1); % reset the random number generator
randprocess=rand(N,T);
iter=0;

for i = 1:N
for t = 2:T

        % Given a number from the uniform distribution...
        rp = randprocess(i,t) ;

        % Depending on the transition matrix, assign individual to a
        % particular state
        iz = 1;
        for xx=1:1:length(prob)
            if rp>sum(prob(sZ(i,t-1),1:xx)) 
                iz = xx+1;
            end
        end

       sZ(i,t) = iz ;
       state_sim(i,t) = theta(iz) ;
       
       Sdec_squeezed = squeeze(kap(Tg(iz,:))) ;
       % there is something wrong here!
       capital_sim(i,t) = interp1q(kap',Sdec_squeezed',capital_sim(i,t-1)) ;

       iteri = iter+i;
       itert = iter+t;
       disp('agent      time')
       disp ([iteri  itert])
end
end


SS=mean(mean(capital_sim));
figure; 
plot(capital_sim(:,:)');
xlim([1 T])
hold on;
plot([1 T], [K K], '-r')
hold on;
plot([1 T], [SS SS], '-b')
title('time series for capital')
ylabel('capital');
xlabel('time');

K- mean(mean(capital_sim))



