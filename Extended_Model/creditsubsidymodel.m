
%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% Model presented in he TBrazilian Slump and the Government-driven 
% Credit Expansion (2020)

%% Initial setups
clear all;
close all;
Parallel=1; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
tic;

%% Toolkit options 

vfoptions.parallel=Parallel;
simoptions.parallel=Parallel;
heteroagentoptions.verbose=1;
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;

%% Parameters

% Preferences 
Params.beta=0.9;% Discount rate

% Firm-level technology
Params.alpha=0.3;  % Capital share
Params.gamma=0.5; % alpha + gamma must be ~= 1
Params.delta=0.05; % Depreciation rate of physical capital
Params.cf=0.1; % Fixed cost of production

% Adjustment cost of capital
Params.adjustcostparam=0.01;

% Entry and Exit
Params.ce=0.4; % Fixed cost of entry 
Params.lambda=0.1; % Probability of firm exit
% lambda is the average observed exit percentage between 2007--2017 
% (https://sidra.ibge.gov.br/Tabela/2718#resultado)
Params.oneminuslambda=1-Params.lambda; % Probability of survival

% Declare discount factors
DiscountFactorParamNames={'beta','oneminuslambda'};
% Declare percentage of entrants
EntryExitParamNames.MassOfNewAgents={'Ne'};
% Exogenous survival probability
EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};

%% States

% The model has three states, one endogenous state (capital), and tow
% exogenous states (productivity and subsidies)

n_s=200;
n_a=100;
% n_psi is two since psi \in {0,1}

%% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=0.02; % Interest rate on earmarked credit
Params.g_ear=0.1; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.

%% Productivity (s)
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

% Earmarked credit grid
psi_grid=[0;1]; % Using this as a {0,1} helps, e.g., add up earmarked capital.
pi_psi=[1,0;0,1];

%% Exogenous states (matrix z)

% Transition matrix 
% Note: considering that productivity and taxes are independent 
n_z=[n_s,length(psi_grid)];
z_grid=[s_grid; psi_grid];

% Transition matrix for the exogenous s and psi variables
pi_z=kron(pi_psi,pi_s);

%% Endogenous state variables

% grid for capital

% steady-state capital without distotions
a_grid = [0 logspace(-2,3,n_a-1)]';

%% Decision variables
%There is no d variable

d_grid=[]; 
n_d=0;


%% Interest rates
% The model has three interest rates
% 1) household
% 2) market
% 3) international

%% 1.Household interest rate

Params.rhhminusdelta=1/Params.beta-1; 
Params.r_hh=Params.rhhminusdelta+Params.delta; 
%% 2.International interest rate
Params.r_international=0.05;

%% 3.Market interest rate
Params.r_market=Params.r_international;

%% Initial guesses and normalization
Params.w=1; % Normalization
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants
%% Potential New Entrants Distribution over the states (k,z)

pistar_s=ones(size(s_grid))/n_s; % Initial guess
for ii=1:10^3 % a while-loop is better practice, but this is easy and I am lazy atm
    pistar_s=(pi_s)'*pistar_s;
end

%% Aspects of the Endogenous entry
% Exit is exogenous with probability lambda
simoptions.agententryandexit=1;
%simoptions.endogenousexit=0;

% Probability of being in the (k, s, psi) category
EntryExitParamNames.DistOfNewAgents={'upsilon'};

% Conditional entry will allow for different productivity cutoffs for 
%new entrants depending on earmarked-vs-nonearmarked.
EntryExitParamNames.CondlEntryDecisions={'ebar'};
% Takes value of one for enter, zero for not-enter. This is just an initial
%guess as the actual decisions are determined as part of general equilibrium.
Params.ebar=ones([n_a,n_z]); 


if Parallel==2
    Params.upsilon = zeros([n_a, n_z],'gpuArray');
else
    Params.upsilon = zeros([n_a, n_z]);
end
Params.upsilon(1,:,:) = kron(pistar_s,[1-Params.g_ear, Params.g_ear]);

% Some checks
disp('upsilon size')
disp(size(Params.upsilon))
disp('sum of upsilon')
disp(sum(Params.upsilon(:)))

%% Return Function
ReturnFn=@(kprime_val, k_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma, cf, adjustcostparam) ExistingFirm_ReturnFn(kprime_val, k_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma, cf, adjustcostparam);
ReturnFnParamNames={'p','w','r_market','r_ear', 'alpha','gamma', 'cf', 'adjustcostparam'}; %It is important that these are in same order as they appear in 'ExistingFirm_ReturnFn'

%% General Equilibrium Equations
%Now define the functions for the General Equilibrium conditions

FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate={};


%Note: length(AggVars) is as for FnsToEvaluate and length(p) is length(n_p)
heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};

% A 'condlentry' general equilibrium condition will take values of greater
% than zero for firms that decide to enter, less than zero for first that
% decide not to enter (or more accurately, after entry decision they draw
% their state, and then decide to cancel/abort their entry).

GEPriceParamNames={'p','Ne'}; 
    
% % Conditional entry condition   (entry - 1st step)
GeneralEqmEqnParamNames(1).Names={'beta'};
GeneralEqmEqn_CondlEntry = @(ValueFn,GEprices,beta) beta*ValueFn-0;

%  Free entry conditions  (entry - 2nd step)
GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce;

GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry};

%% Find equilibrium prices

heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input 

[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1...
    (n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, [],...
    EntryExitParamNames);

Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;

%% Value Function, Policy and Firm Distribution in GE
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames,ReturnFnParamNames);
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,...
    Params, EntryExitParamNames);

%% Solve the partial equilibrium problem of capital market clearance.
% Calculate K_nfa

%Params.K_nfa is the net foreign assets
% That it the total physical capital held by the country of foreign assets.
% Note that an increase in this would represent positive net exports 

% Interest rate for firms is r_market == international interest rate
% Interest rate for households is r_HH, 

% To solve the partial equilibrium:
% STEP 1 - Let K_total=K_firm(r_market) 
% STEP 2 - Let K_HH=K_firm(r_HH)  
% STEP 3 - Now, set K_nfa=K_hh-K_total  


% STEP 1 - Let K_total=K_firm(r_market) 
Params.r_market=Params.r_market;
FnsToEvaluateParamNames(1).Names={};
% Find Aggregate Capital with r = r_market
FnsToEvaluateFn_Kfirm = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass)...
    aprime_val; 
FnsToEvaluate={FnsToEvaluateFn_Kfirm};
K_total=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, [],simoptions,EntryExitParamNames);

% STEP 2 - Let K_HH=K_firm(r_HH)  
% Params.rhhminusdelta=1/Params.beta-1, that is general eqm result in 
% complete market models
% Params.r_hh=Params.rhhminusdelta+Params.delta;
Params.r_market=Params.r_hh;
% Policy function with r = r_HH
[~,Policy_rhh]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames,ReturnFnParamNames); 
% Find Aggregate Capital with r = r_HH
K_hh=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy_rhh,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, [],simoptions,EntryExitParamNames);


%% 
%Set r_market back to it's actual value
Params.r_market=Params.r_international;

% 2C: Calculate K_nfe
K_nfa=K_hh-K_total;

%%

%
fprintf('\n')
fprintf('Capital Market Outcomes: \n')
fprintf('The interest rate for household is %8.2f, the international interest rate is r%8.2f, and the interest rate on earmarked credit is %8.2f \n', 100*Params.r_hh, 100*Params.r_international, 100*Params.r_ear)
fprintf('Households hold capital of K_hh=%8.4f, while firms use capital of K_total=%8.4f; the difference is made up by net-foreign capital of K_nfa=%8.4f \n',K_hh, K_total, K_nfa)
fprintf('\n')

toc;
