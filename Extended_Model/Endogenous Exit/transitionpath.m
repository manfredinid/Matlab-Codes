%%

T=50 % number of time periods to transition path

Params=Params_initial;

ParamPath=Params.taurate_final*ones(T,1);
ParamPathNames={'tau'};

transpathoptions.verbose=1
transpathoptions.GEnewprice=2


%% Transition Path
transpath_shootingalgo=0;
simoptions.agententryandexit=1;
simoptions.endogenousexit=1;
transpathoptions.agententryandexit=1;
transpathoptions.parallel=2;
transpathoptions.exoticpreferences=0;
transpathoptions.weightscheme=1;


%% Return Function
DiscountFactorParamNames={'beta'};


ReturnFn=@(kprime_val, k_val,s_val, psi_val, p,w,r_market,r_ear,...
    alpha,gamma,delta, cf, adjustcostparam) ExistingFirm_ReturnFn(kprime_val,...
    k_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma,delta, cf, adjustcostparam);
ReturnFnParamNames={'p','w','r_market','r_ear', 'alpha','gamma','delta',...
    'cf', 'adjustcostparam'}; %It is important that these are in same order as they appear in 'ExistingFirm_ReturnFn'



%%
FnsToEvaluateParamNames(2).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma);  

FnsToEvaluate={FnsToEvaluateFn_output}

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

AggVars

% We need to give an initial guess for the price path on interest rates
PricePath0_p=[linspace(Params_initial.p, Params_final.p, floor(T/2))'; Params_final.p*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0_Ne=[linspace(Params_initial.Ne, Params_final.Ne, floor(T/2))'; Params_final.Ne*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0_r=[linspace(Params_initial.r, Params_final.r, floor(T/2))'; Params_final.r*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0=[PricePath0_p, PricePath0_Ne]; % PricePath0 is matrix of size T-by-'number of prices'
PricePathNames={'p','Ne'};
%%
% Rewrite the General Eqm conditions as rules for updating the price
heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};

 % Alternative attempt, based on minimizing weighted sum of squares.
    transpathoptions.weightsforpath=ones(T,2); % Same size as PricePath

GEPriceParamNames={'ebar'}; 
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate={};
heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};

GeneralEqmEqnParamNames(1).Names={'beta'};
GeneralEqmEqn_CondlEntry = @(ValueFn,GEprices,beta) beta*ValueFn-0;
GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.


 heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};
 GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry}; 
%%
[PricePath]= TransitionPath_Case1(PricePath0, PricePathNames,...
 ParamPath, ParamPathNames, T, V_final, StationaryDist_initial,...
   n_d,n_a, n_z, pi_z,d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate,...
   GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    transpathoptions, vfoptions, simoptions, EntryExitParamNames);

    save ./SavedOutput/PricePath