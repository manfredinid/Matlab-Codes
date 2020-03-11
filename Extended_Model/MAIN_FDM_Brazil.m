
%% Stationary Equilibrium and Transition path for the Firm Dynamics Model

%clear all;
%close all;

Parallel=0 % 2 for GPU, 1 for parallel CPU, 0 for single CPU.
SkipInitialFinal= 1 % 1 to SKIP transition path

%% Endogenous and Exogenous States
n_s= 5; % number of firm-specific Productivity level
n_psi = 10; % number of credit tax 
n_a=20; % grid size for capital

%% Stacionary Equilibrium
if SkipInitialFinal==1
fprintf(2,'\nStacionary Equilibrium\n')
    
%Policy parameters
Params.gcost=0.01;   
psi_grid = linspace(-10,10,n_psi)';

% Initial guesses
Params.p=1; % output pricecap
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Parameters_FDM_Brazil;
% Stationary Equilibrium and Results
SS_FDM_Brazil;


%% Transition Path
else
   fprintf(2,'\nTransition Path\n') 

%% Transition

% credit tax
Params.psi_initial = linspace(-1,1,n_psi)';
Params.psi_final = linspace(0,5,n_psi)';

% adjustment cost parameter
Params.gcost_initial=0.01;
Params.gcost_final = 0.05;
%% Initial Period  

%Policy parameters
Params.gcost=Params.gcost_initial;   
psi_grid = Params.psi_initial;

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

disp('Calculating various equilibrium objects')
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

%Policy parameters
Params.gcost=Params.gcost_final;   
psi_grid = Params.psi_final;


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

disp('Calculating various equilibrium objects')
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
    [V_final,Policy_final]=ValueFnIter_Case1(V0, n_d,n_a,n_z,[],a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);

StationaryDist_final=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions, Params, EntryExitParamNames);

Params_final=Params;
save ./SavedOutput/HopenhaynRogerson1993_final.mat Params_final...
    V_final Policy_final StationaryDist_final

%% General Equilibrium Transition Path

T=50 % number of time periods to transtion path

Params=Params_initial;


transpathoptions.parallel=0;
transpath_shootingalgo=0;
vfoptions.endogenousexit=0;
transpathoptions.agentexit=0;
transpathoptions.agententry=1;
%%
ParamPath=Params.gcost_final*ones(T,1);
ParamPathNames={'gcost'};

% We need to give an initial guess for the price path on interest rates
PricePath0_p=[linspace(Params_initial.p, Params_final.p, floor(T/2))'; Params_final.p*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0_Ne=[linspace(Params_initial.Ne, Params_final.Ne, floor(T/2))'; Params_final.Ne*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0=[PricePath0_p, PricePath0_Ne]; % PricePath0 is matrix of size T-by-'number of prices'
PricePathNames={'p','Ne'};

% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.specialgeneqmcondn={0,'entry'};

 % Alternative attempt, based on minimizing weighted sum of squares.
    transpathoptions.GEnewprice=2;
    transpathoptions.weightsforpath=ones(T,2); % Same size as PricePath

GeneralEqmEqnParamNames(1).Names={};
GeneralEqmEqn_GoodsMarket = @(AggVars,GEprices) 1-AggVars;

GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.

GeneralEqmEqns={GeneralEqmEqn_GoodsMarket,GeneralEqmEqn_Entry};   

transpathoptions.weightscheme=1
transpathoptions.verbose=1


[PricePath]=TransitionPath_Case1(PricePath0, PricePathNames, ParamPath,....
    ParamPathNames, T, V_final, StationaryDist_initial, n_d, n_a, n_z, pi_z,...
    d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params,...
    DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames,...
    GeneralEqmEqnParamNames,transpathoptions, vfoptions, simoptions, EntryExitParamNames);

end