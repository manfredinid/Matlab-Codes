%% Transition Path
clear FnsToEvaluateParamNames


transpath_shootingalgo=0;
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;
transpathoptions.agententryandexit=1;
transpathoptions.parallel=2;
transpathoptions.exoticpreferences=0;
transpathoptions.weightscheme=1;


%% Return Function
DiscountFactorParamNames={'beta'};

ReturnFn=@(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,delta,taurate,subsidyrate, cf, gcost)...
capital_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma, delta, taurate,subsidyrate, cf, gcost);

ReturnFnParamNames={ 'p','r', 'alpha','gamma', 'delta','taurate','subsidyrate', 'cf', 'gcost'};

%%

FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_nbar =@(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate)...
((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));

FnsToEvaluate={FnsToEvaluateFn_nbar};

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel, simoptions, EntryExitParamNames);

AggVars

% We need to give an initial guess for the price path on interest rates
PricePath0.p=[linspace(Params_initial.p, Params_final.p, floor(T/2))'; Params_final.p*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
PricePath0.Ne=[linspace(Params_initial.Ne, Params_final.Ne, floor(T/2))'; Params_final.Ne*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePath0.r=[linspace(Params_initial.r, Params_final.r, floor(T/2))'; Params_final.r*ones(T-floor(T/2),1)]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePath0=[PricePath0_p, PricePath0_Ne]; % PricePath0 is matrix of size T-by-'number of prices'
% PricePathNames={'p','Ne'};

%%
% Rewrite the General Eqm conditions as rules for updating the price
transpathoptions.specialgeneqmcondn={'entry',0};

 % Alternative attempt, based on minimizing weighted sum of squares.
 transpathoptions.weightsforpath=ones(T,2); % Same size as PricePath

GeneralEqmEqnParamNames(1).Names={'beta','ce'};
GeneralEqmEqn_Entry2 = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce;
GeneralEqmEqnParamNames(2).Names={};
GeneralEqmEqn_GoodsMarket2 = @(AggVars,GEprices) 1-AggVars;

GeneralEqmEqns={GeneralEqmEqn_Entry2,GeneralEqmEqn_GoodsMarket2};   

%%
[PricePath]= TransitionPath_Case1(PricePath0, ParamPath, T, V_final, StationaryDist_initial,...
   n_d,n_a, n_z, pi_z,d_grid,a_grid,z_grid, ReturnFn, FnsToEvaluate,...
   GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    transpathoptions, vfoptions, simoptions, EntryExitParamNames);


