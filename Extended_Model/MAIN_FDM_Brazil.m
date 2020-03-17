
%% Stationary Equilibrium and Transition path for the Firm Dynamics Model

clear all;
close all;

Parallel=1 % 2 for GPU, 1 for parallel CPU, 0 for single CPU.
SkipInitialFinal= 1 % 1 to SKIP transition path


tic;
%% Endogenous and Exogenous States
n_s= 10; % number of firm-specific Productivity level
n_psi = 3; % number of credit tax 
n_a=30; % grid size for capital

%% Stacionary Equilibrium

if SkipInitialFinal==1
fprintf(2,'\nStacionary Equilibrium\n')
    
%Policy parameters
Params.gcost=0.01;   
% Distortions
Params.taurate=0.1; % This is the rate for the tax.
Params.subsidyrate=0.01; % This is the rate for the subsidy.

% subsidy-tax distribution (new entrants)
%psi_grid = linspace(-1,1,n_psi)'; % Incumbest first draws
%psi_dist =betarnd(.5,.4, 1, n_psi); % Entrants probability distribution
psi_grid = [-1; 0; 1]; % Incumbest first draws
psi_dist = [0.5; 0; 0.5]; % Entrants probability distribution
% Why I did not use just 3 values in the psi_grid?
% Because this way I have more control over the probability distribution of psi
% Maybe a 3 variables grid was best - it was not

% Initial guesses
Params.p=1; % output pricecap
Params.Ne=0.5; % total mass of new entrants
%%
% Parameters and initialization options
Parameters_FDM_Brazil;
% Stationary Equilibrium and Results
SS_FDM_Brazil;


%% Transition Path
else
   fprintf(2,'\nTransition Path\n') 



% Distortions

% INITIAL
Params.taurate_initial=0.2; % This is the rate for the tax.
Params.subsidyrate_initial=0.2; % This is the rate for the subsidy.
Params.gcost_initial=0.01;

psi_grid_initial = [-1; 0; 1]; % Incumbest first draws
psi_dist_initial =[0; 0; 1]; 

% FINAL
Params.taurate_final=0.6; % This is the rate for the tax.
Params.subsidyrate_final=0.2; % This is the rate for the subsidy.
Params.gcost_final = 0.01;

% psi_grid_final should be the same in the initial and final period
% psi_dist_final should change based on the policy
psi_grid_final = [-1; 0; 1]; % Incumbest first draws
psi_dist_final = [0; 0; 1]; 

% Why I did not use just 3 values in the psi_grid?

%% Initial Period  

  fprintf('\nInitial Period\n') 
  
%Policy parameters
Params.gcost=Params.gcost_initial;   
Params.subsidyrate=Params.subsidyrate_initial;
Params.taurate=Params.taurate_initial;

psi_grid = psi_grid_initial;
psi_dist = psi_dist_initial;

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Parameters_FDM_Brazil;

% Find equilibrium prices
heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,p_eqm_index,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(V0,...
    n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
Params_initial=Params;


% Value Function, Policy and Firm Distribution in GE

disp('Calculating value function and policy function')
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
  [V_initial,Policy_initial]=ValueFnIter_Case1(V0, n_d,n_a,n_z,[],a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);

StationaryDist_initial=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions, Params, EntryExitParamNames);

Params_initial=Params;

    save ./SavedOutput/TPDynamics_initial.mat...
        Params_initial V_initial Policy_initial StationaryDist_initial

%% Final Period

  fprintf('\nFinal Period\n') 

%Policy parameters
Params.gcost=Params.gcost_final;   
Params.subsidyrate=Params.subsidyrate_final;
Params.taurate=Params.taurate_final;

psi_grid = psi_grid_initial;
psi_dist = psi_dist_final;

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Parameters_FDM_Brazil;

% Find equilibrium prices
heteroagentoptions.verbose=1;

n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,p_eqm_index,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(V0,...
    n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);


% Value Function, Policy and Firm Distribution in GE

disp('Calculating value function and policy function')
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
    [V_final,Policy_final]=ValueFnIter_Case1(V0, n_d,n_a,n_z,[],a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);

StationaryDist_final=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions, Params, EntryExitParamNames);

Params_final=Params;

 save ./SavedOutput/TPDynamics_final.mat...
        Params_final V_final Policy_final StationaryDist_final

%% General Equilibrium Transition Path
 
 fprintf('\nGeneral Equilibrium Transition Path\n') 
 
T=50 % number of time periods to transition path

Params=Params_initial;


transpathoptions.parallel=2;
transpath_shootingalgo=0;
transpathoptions.agententry=1;

%% Return Function
ReturnFn=@(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,delta,taurate,subsidyrate, cf, gcost)...
FDM_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,delta,taurate,subsidyrate, cf, gcost);

ReturnFnParamNames={ 'p','r', 'alpha','gamma', 'delta','taurate','subsidyrate', 'cf', 'gcost'};
%%
ParamPath=Params.taurate_final*ones(T,1);
ParamPathNames={'tau'};

% We need to give an initial guess for the price path on interest rates
PricePath0_p=[linspace(Params_initial.p, Params_final.p, floor(T/2))'; Params_final.p*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0_Ne=[linspace(Params_initial.Ne, Params_final.Ne, floor(T/2))'; Params_final.Ne*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0=[PricePath0_p, PricePath0_Ne]; % PricePath0 is matrix of size T-by-'number of prices'
PricePathNames={'p','Ne'};

% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.specialgeneqmcondn={'entry',0};

 % Alternative attempt, based on minimizing weighted sum of squares.
    transpathoptions.GEnewprice=2;
    transpathoptions.weightsforpath=ones(T,2); % Same size as PricePath

GeneralEqmEqnParamNames(1).Names={'beta','ce'};
GeneralEqmEqn_Entry2 = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce;
GeneralEqmEqnParamNames(2).Names={};
GeneralEqmEqn_GoodsMarket2 = @(AggVars,GEprices) 1-AggVars;



GeneralEqmEqns={GeneralEqmEqn_GoodsMarket2,GeneralEqmEqn_Entry2};   

transpathoptions.weightscheme=1
transpathoptions.verbose=1


[PricePath]=TransitionPath_Case1(PricePath0, PricePathNames,...
 ParamPath, ParamPathNames, T, V_final, StationaryDist_initial,...
  n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate,...
   GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    transpathoptions, vfoptions, simoptions, EntryExitParamNames);

end

toc;