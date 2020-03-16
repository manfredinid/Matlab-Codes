%% Credit Model with Firm Dynamics


%clear all;
%close all;
% Second version

%Parallel=0 % 2 for GPU, 1 for parallel CPU, 0 for single CPU.

rng('default') % For reproducibility

%% Toolkit options
tauchenoptions.parallel=Parallel;

mcmomentsoptions.T=10^4;
mcmomentsoptions.Tolerance=10^(-9);
mcmomentsoptions.parallel=tauchenoptions.parallel;

vfoptions.parallel=Parallel;

simoptions.burnin=10^4;
simoptions.simperiods=10^10; % if iterate=0 then simperiod=10^6 
simoptions.iterate=1;
simoptions.parallel=Parallel; 

simoptions.maxit=10^4;

heteroagentoptions.verbose=1;

%% Parameters Calibration

% Preferences 
Params.beta=0.9;% Discount rate

% Firm-level technology
Params.alpha=0.3;  % Capital share
Params.gamma=0.5; % alpha + gamma must be ~= 1
Params.delta=0.05; % Depreciation rate of physical capital
Params.cf=1; % Fixed cost of production

% Entry and Exit
Params.ce=1; % Fixed cost of entry 
Params.lambda=0.2; % Probability of firm exit
% lambda is the average observed exit percentage between 2007--2017 
% (https://sidra.ibge.gov.br/Tabela/2718#resultado)
Params.oneminuslambda=1-Params.lambda; % Probability of survival

% Distortions
%Params.taurate=5; % This is the rate for the tax.
%Params.subsidyrate=5; % This is the rate for the subsidy.
%Params.gcost=0.01; % capital adjustment cost parameter

% Initial guesses
%Params.p=1; % output price
%Params.Ne=0.5; % total mass of new entrants

% Declare discount factors
DiscountFactorParamNames={'beta','oneminuslambda'};
% Declare percentage of entrants
EntryExitParamNames.MassOfNewAgents={'Ne'};
% Exogenous survival probability
EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};

%% Steady-state interest rate

Params.i=1/Params.beta-1; % gross capital return
Params.r=Params.i+Params.delta; % net capital return

%% Exogenous state variables

%n_s= 20; % firm-specific Productivity level
%n_psi = 5; % credit tax 


% Exogenous AR(1) process on (log) productivity
% logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)
Params.rho=0.93; 
Params.sigma_logz=sqrt(0.53); 
Params.sigma_epsilon=sqrt((1-Params.rho)*((Params.sigma_logz)^2));
Params.a=0.098; 

tauchenoptions.parallel=Parallel;
Params.q=2; % Hopenhayn & Rogerson (1993) do not report (based on Table 4 is seems something around q=4 is used, otherwise don't get values of z anywhere near as high as 27.3. (HR1993 have typo and call the column 'log(s)' when it should be 's') 
[s_grid, pi_s]=TauchenMethod(Params.a,Params.sigma_epsilon^2,Params.rho,n_s,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
s_grid=exp(s_grid);

% Tax credit
%psi_grid = linspace(-1,1,n_psi)';


% Transition matrix 
% Note: considering that productivity and taxes are independent 
n_z=[n_s,length(psi_grid)];
z_grid=[s_grid; psi_grid];


% transition matrix for the exogenous z and psi variables
pi_z=kron( pi_s,eye(prod(n_psi)))';

% Check transition matrix
%for ii = 1: length(pi_z)
%A = round(sum(pi_z(:,ii)),5);
%if A == 1
%else
%   error('transition matrix sum is not one')
%end
%end
pi_z=pi_z';

%% Endogenous state variables

% grid for capital
%n_a=50;

% steady-state capital without distotions
%%%%% The grid is like the one in the Aiygari example
k_ss = (Params.alpha/Params.r)^(1-Params.gamma/1-Params.gamma-Params.alpha)*...
    (Params.gamma)^(Params.gamma/1-Params.alpha-Params.gamma);
nk1 = floor(n_a/3); nk2=floor(n_a/3); nk3=n_a-nk1-nk2;
a_grid = sort([linspace(0,k_ss,nk1),linspace(k_ss+0.0001,3*k_ss,nk2),...
      linspace(3*k_ss+0.0001,15*k_ss,nk3)])';

%% Decision variables
%There is no d variable

d_grid=[]; 
n_d=0;

%% Check endogenous, exogenous and decision variables

disp('sizes')
disp('vector(s) of endogenous state variables')
disp(n_a)
disp('vector(s) of exogenous state variable')
disp(n_z)
disp('vector(s) of decision variabes')
disp(n_d)

%% Potential New Entrants Distribution over the states (s, psi, k)

% productivity (exogenous state)

logn = lognrnd(1,0.5,1,n_s);
cumsum_pistar_s = cumsum(logn./sum(logn));
pistar_s=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)]);

% credit tax (exogenous state)
%beta = betarnd(.5,.4, 1, n_psi);
cumsum_pistar_psi = cumsum(psi_dist./sum(psi_dist))';
pistar_psi =(cumsum_pistar_psi-[0,cumsum_pistar_psi(1:end-1)])';

% capital (endogenous state)
pistar_k = [zeros(n_a-1,1);1];
cumsum_pistar_k = cumsum(pistar_k);

if (abs(1-round(sum(pistar_psi),2)) || abs(1-round(sum(pistar_psi),2))||abs(1-sum(pistar_k)) > 1e-5)
   error('Draws are NOT a PMD.')
end


figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-|--|:','DefaultLineLineWidth',1);
subplot(3,1,1);
plot(psi_grid,cumsum_pistar_psi,'r')
hold on;
line([0,0], [0 1])
title('Potential draws for psi')
subplot(3,1,2);
plot(s_grid,cumsum_pistar_s,'r')
title('Potential draws for s')
subplot(3,1,3);
plot(a_grid,cumsum_pistar_k,'r')
title('Potential draws for k')

%% Return Function
ReturnFn=@(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,delta,taurate,subsidyrate, cf, gcost)...
FDM_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,delta,taurate,subsidyrate, cf, gcost);

ReturnFnParamNames={ 'p','r', 'alpha','gamma', 'delta','taurate','subsidyrate', 'cf', 'gcost'};

%% CHECK (to be erase)
if vfoptions.parallel==2
    V0=zeros([n_a,n_z],'gpuArray');
else
    V0=zeros([n_a,n_z]);
end
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,...
    a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, vfoptions);

figure;
surf(squeeze(V(:,1,:)))
 
figure;
surf(squeeze(Policy(1,:,1,:)))

%% Aspects of the Endogenous entry
% Exit is exogenous with probability lambda
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;

% Probability of being in the (s, psi) category
EntryExitParamNames.DistOfNewAgents={'upsilon'};

pistar_psi_s=pistar_s'.*(pistar_psi)';
Params.upsilon=NaN(n_a,n_s,n_psi);
 for n=1:n_a
    Params.upsilon(n,:,:)=pistar_psi_s.*pistar_k(n);
 end

 
disp('upsilon size')
disp(size(Params.upsilon))

disp('sum of upsilon')
disp(sum(Params.upsilon(:)))

%% CHECK (to be erased)
simoptions.parallel=Parallel

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions,Params,EntryExitParamNames);

figure;
surf(squeeze(StationaryDist.pdf(:,1,:)))



%%


%Use the toolkit to find the equilibrium price index
GEPriceParamNames={'p'};

%FnsToEvaluateParamNames(1).Names={};
%FnsToEvaluate={};

heteroagentoptions.specialgeneqmcondn={0,'entry'};

FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_nbar =@(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate)...
(((1-taurate*z2_val)*p*z1_val*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)); 
FnsToEvaluate={FnsToEvaluateFn_nbar};
%%
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

AggVars

%%
GEPriceParamNames={'p', 'Ne'}; 
GeneralEqmEqnParamNames(1).Names={};
GeneralEqmEqn_GoodsMarket = @(AggVars,GEprices) 1-AggVars;

GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.

GeneralEqmEqns={GeneralEqmEqn_GoodsMarket,GeneralEqmEqn_Entry};

%%

