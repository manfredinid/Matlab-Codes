%% Credit Model with Firm Dynamics

rng('default') % For reproducibility

%% Toolkit options
tauchenoptions.parallel=Parallel;

mcmomentsoptions.T=10^4;
mcmomentsoptions.Tolerance=10^(-9);
mcmomentsoptions.parallel=tauchenoptions.parallel;

vfoptions.parallel=Parallel;

simoptions.burnin=10^4;
simoptions.simperiods=10^6; % if iterate=0 then simperiod=10^6 
simoptions.iterate=0;
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
Params.cf=0.1; % Fixed cost of production

% Entry and Exit
Params.ce=0.4; % Fixed cost of entry 
Params.lambda=0.1; % Probability of firm exit
% lambda is the average observed exit percentage between 2007--2017 
% (https://sidra.ibge.gov.br/Tabela/2718#resultado)
Params.oneminuslambda=1-Params.lambda; % Probability of survival

% Distortions
% Select in the main file

% Initial guesses
% Select in the main file

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

% Exogenous AR(1) process on (log) productivity
% logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)
Params.rho=0.93; 
Params.sigma_logz=sqrt(0.53); 
Params.sigma_epsilon=sqrt((1-Params.rho)*((Params.sigma_logz)^2));
Params.a=0.098; 

tauchenoptions.parallel=Parallel;
Params.q=2; 
[s_grid, pi_s]=TauchenMethod(Params.a,Params.sigma_epsilon^2,Params.rho,n_s,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
s_grid=exp(s_grid);

% Tax credit
% Select in the main file


% Transition matrix 
% Note: considering that productivity and taxes are independent 
n_z=[n_s,length(psi_grid)];
z_grid=[s_grid; psi_grid];


% transition matrix for the exogenous z and psi variables
pi_z=kron( eye(prod(n_psi)),pi_s);


%% Endogenous state variables

% grid for capital
% Select in the main file

% steady-state capital without distotions
  a_grid = [0 logspace(-2,3,n_a-1)]';

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
%logn = lognrnd(1,0.5,1,n_s);
%cumsum_pistar_s = cumsum(logn./sum(logn));
%pistar_s=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)]);
pistar_s=ones(size(s_grid))/n_s; % Initial guess
dist=1;
while dist>10^(-9)
    pistar_s_old=pistar_s;
    pistar_s=(pi_s')*pistar_s;
    dist=max(abs(pistar_s-pistar_s_old));
end


% credit tax (exogenous state)
pistar_psi =psi_dist;

if (abs(1-round(sum(pistar_psi),2)) ||abs(1-sum(pistar_s)) > 1e-5)
   error('Draws are NOT a PMD.')
end


figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-|--|:','DefaultLineLineWidth',1);
subplot(2,1,1);
plot(psi_grid,pistar_psi,'r')
hold on;
line([0,0], [0 1])
title('Potential draws for psi')
subplot(2,1,2);
plot(s_grid,pistar_s,'r')
title('Potential draws for s')


%% Return Function
ReturnFn=@(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,delta,taurate,subsidyrate, cf, gcost)...
capital_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma, delta, taurate,subsidyrate, cf, gcost);

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

figure; surf(squeeze(Policy(1,:,:,2)))
xlabel('productivity')
ylabel('capital')
title('For psi=0')

%% Aspects of the Endogenous entry
% Exit is exogenous with probability lambda
simoptions.agententryandexit=1;
%simoptions.endogenousexit=0;

% Probability of being in the (s, psi) category
EntryExitParamNames.DistOfNewAgents={'upsilon'};

Params.upsilon = zeros([n_a, n_z],'gpuArray');
Params.upsilon(1,:,:) = kron(pistar_s,(pistar_psi)');


disp('upsilon size')
disp(size(Params.upsilon))

disp('sum of upsilon')
disp(sum(Params.upsilon(:)))

%% CHECK (to be erased)
simoptions.parallel=Parallel

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions,Params,EntryExitParamNames);

StationaryDist.mass

% only k=0 has positive probabilities
figure; surf(squeeze(StationaryDist.pdf(:,:,2)))

%% Use the toolkit to find the equilibrium price index
GEPriceParamNames={'p'};

FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate={};

heteroagentoptions.specialgeneqmcondn={0,'entry'};

FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_nbar =@(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate)...
((z1_val*p*gamma))^(1/(1-gamma))*aprime_val^(alpha/(1-gamma));




FnsToEvaluate={FnsToEvaluateFn_nbar};
%%
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

AggVars
Output.perK=AggVars/StationaryDist.mass;
Output.perK

%%
GEPriceParamNames={'p', 'Ne'}; 
GeneralEqmEqnParamNames(1).Names={};
GeneralEqmEqn_nbar = @(AggVars,p) 1-AggVars;

GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,p,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.

GeneralEqmEqns={GeneralEqmEqn_nbar,GeneralEqmEqn_Entry};

%%

