
%% Stationary Equilibrium and Transition path for the Firm Dynamics Model

clear all;
close all;

Parallel=0 % 2 for GPU, 1 for parallel CPU, 0 for single CPU.
SkipInitialFinal= 0 % 1 to SKIP transition path

%% Endogenous and Exogenous States
n_s= 20; % number of firm-specific Productivity level
n_psi = 5; % number of credit tax 
n_a=50; % grid size for capital

%% Stacionary Equilibrium
if SkipInitialFinal==1
fprintf(2,'\nStacionary Equilibrium\n')
    
%Policy parameters
Params.gcost=0.01;   
psi_grid = linspace(-1,1,n_psi)';

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Parameters_dynamics_credit;
% Stationary Equilibrium and Results
SS_dynamics_credit;


%% Transition Path
else
   fprintf(2,'\nTransition Path\n') 
%% Initial Period  

%Policy parameters
Params.gcost=0.01;   
psi_grid = linspace(-1,1,n_psi)';

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Parameters_dynamics_credit;

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
    save ./SavedOutput/TPDynamics_initial.mat...
        Params_initial V_initial Policy_initial ExitPolicy_initial StationaryDist_initial

%% Final Period

%Policy parameters
Params.gcost=0.01;   
psi_grid = linspace(0,5,n_psi)';

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Parameters_dynamics_credit;

% Find equilibrium prices
heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,p_eqm_index,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(V0,...
    n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
save ./SavedOutput/TPDynamics_final.mat...
        Params_initial V_initial Policy_initial ExitPolicy_initial StationaryDist_initial

%% General Equilibrium Transition Path

T=50 % number of time periods to transtion path

Params=Params_initial;

transpathoptions.agententryandexit=1
%%%%%%%%% ERROR transition path just work with parallel = 2
%transpathoptions.parallel=0;
transpath_shootingalgo=0

%%
ParamPath=Params.psi_final*ones(T,1);
ParamPathNames={'psi'};

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
    
 GEPriceParamNames={'p', 'Ne'}; 
GeneralEqmEqnParamNames(1).Names={};
GeneralEqmEqn_LabourMarket = @(AggVars,GEprices) 1-AggVars;

GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.

GeneralEqmEqns={GeneralEqmEqn_LabourMarket,GeneralEqmEqn_Entry};   

transpathoptions.weightscheme=1
transpathoptions.verbose=1
[PricePath]=TransitionPath_Case1(PricePath0, PricePathNames, ParamPath,....
    ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z,...
    d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params,...
    DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames,...
    GeneralEqmEqnParamNames,transpathoptions)%, vfoptions%, simoptions, EntryExitParamNames);

end