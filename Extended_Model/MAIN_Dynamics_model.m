
%% Stationary Equilibrium and Transition path for the Firm Dynamics Model

clear all;
close all;

Parallel=0 % 2 for GPU, 1 for parallel CPU, 0 for single CPU.
SkipInitialFinal= 1 % 1 to SKIP transition path

%% Endogenous and Exogenous States
n_s= 20; % firm-specific Productivity level
n_psi = 5; % credit tax 
n_a=50; % grid for capital

%% Stacionary Equilibrium
if SkipInitialFinal==1
 
%Policy parameters
Params.gcost=0.01;   
psi_grid = linspace(-1,1,n_psi)';

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Dynamics_credit_model;
% Stationary Equilibrium and Results
stationary;


%% Transition Path
else
%% Initial Period  

%Policy parameters
Params.gcost=0.01;   
psi_grid = linspace(-1,1,n_psi)';

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Dynamics_credit_model

% Find equilibrium prices
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
psi_grid = linspace(-1,1,n_psi)';

% Initial guesses
Params.p=1; % output price
Params.Ne=0.5; % total mass of new entrants

% Parameters and initialization options
Dynamics_credit_model

% Find equilibrium prices
disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,p_eqm_index,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(V0,...
    n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
save ./SavedOutput/TPDynamics_final.mat...
        Params_initial V_initial Policy_initial ExitPolicy_initial StationaryDist_initial
end