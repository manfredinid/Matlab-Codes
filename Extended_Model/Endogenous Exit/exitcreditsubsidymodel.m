
%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION (endogenous exit)
% Model presented in the Brazilian Slump and the Government-driven 
% Credit Expansion (2020)

%% Initial setups
%clear all;
%close all;
%Parallel=1; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
%tic;
%rng(1);

%% Toolkit options 

vfoptions.parallel=Parallel;
simoptions.parallel=Parallel;
heteroagentoptions.verbose=1;
simoptions.agententryandexit=1;
simoptions.endogenousexit=1;


%% Parameters

% Preferences 
Params.beta=0.9798;% Discount rate

% Firm-level technology
Params.alpha=0.399;  % Capital share
Params.gamma=0.491; % alpha + gamma must be ~= 1
Params.delta=0.025; % Depreciation rate of physical capital
Params.cf=5; % Fixed cost of production
%3

Params.w=1; % Normalization
Params.p=0.3549; % output price
% Adjustment cost of capital
Params.adjustcostparam = 3.219;

% Entry and Exit
Params.ce=5*Params.cf; % Fixed cost of entry 
%4
%% States

% The model has three states, one endogenous state (capital), and tow
% exogenous states (productivity and subsidies)

n_s=15;
n_a=350;
% n_psi is two since psi \in {0,1}

%% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

%Params.r_ear=0.02; % Interest rate on earmarked credit
%Params.g_ear=0.4; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.

%% Productivity (s)
% Exogenous AR(1) process on (log) productivity
% logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)
Params.rho=0.9; 
Params.sigma_logz=0.13; 
Params.sigma_epsilon=sqrt((1-Params.rho)*((Params.sigma_logz)^2));
Params.a=0.04; 

tauchenoptions.parallel=Parallel;
Params.q=5; 
%  q is max number of std devs from mean
% max in s_grid has to be log(2.75)=1.0116
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
a_grid = [0 logspace(0.0001,6.28,n_a-1)]'; 


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
%Params.r_international=0.15;

%% 3.Market interest rate
Params.r_market=Params.r_international;



%% Potential New Entrants Distribution over the states (k,z)

pistar_s=ones(size(s_grid))/n_s; % Initial guess
dist=1;
while dist>10^(-7)
    pistar_s_old=pistar_s;
    pistar_s=(pi_s)'*pistar_s;
    dist=max(abs(pistar_s-pistar_s_old));
end

%% Aspects of the Endogenous entry
% Exit is exogenous with probability lambda
DiscountFactorParamNames={'beta'};
vfoptions.endogenousexit=1;



%% Return Function
ReturnFn=@(kprime_val, k_val,s_val, psi_val, p,w,r_market,r_ear,...
    alpha,gamma,delta, cf, adjustcostparam) ExistingFirm_ReturnFn(kprime_val,...
    k_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma,delta, cf, adjustcostparam);
ReturnFnParamNames={'p','w','r_market','r_ear', 'alpha','gamma','delta',...
    'cf', 'adjustcostparam'}; %It is important that these are in same order as they appear in 'ExistingFirm_ReturnFn'

vfoptions.ReturnToExitFn=@(k_val,s_val, psi_val)0; 
vfoptions.ReturnToExitFnParamNames={};

[V,Policy,ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid,...
    pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);


%%

% Probability of being in the (k, s, psi) category
EntryExitParamNames.DistOfNewAgents={'upsilon'};

EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};
Params.oneminuslambda=1-ExitPolicy;

% Conditional entry will allow for different productivity cutoffs for 
%new entrants depending on earmarked-vs-nonearmarked.
EntryExitParamNames.CondlEntryDecisions={'ebar'};
% Takes value of one for enter, zero for not-enter. This is just an initial
%guess as the actual decisions are determined as part of general equilibrium.
Params.ebar=zeros([n_a,n_z]); 

EntryExitParamNames.MassOfNewAgents={'Ne'};

if Parallel==2
    Params.upsilon = zeros([n_a, n_z],'gpuArray');
else
    Params.upsilon = zeros([n_a, n_z]);
end
Params.upsilon(1,:,:) = kron(pistar_s,[1-Params.g_ear, Params.g_ear]);

% Some checks
disp('upsilon size')
disp(size(Params.upsilon))
disp('nansum of upsilon')
disp(nansum(Params.upsilon(:)))

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,...
    Params, EntryExitParamNames);
%% General Equilibrium Equations
%Now define the functions for the General Equilibrium conditions
GEPriceParamNames={'ebar'}; 
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluate={};
heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};
%FnsToEvaluateParamNames(1).Names={'p', 'w','alpha','gamma'};
%FnsToEvaluateFn_1 = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
%    z1_val*(aprime_val^alpha)*...
%    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma); 
%FnsToEvaluate={FnsToEvaluateFn_1};

%simoptions.keeppolicyonexit=1;
%AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
%    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
%    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

% A 'condlentry' general equilibrium condition will take values of greater
% than zero for firms that decide to enter, less than zero for first that
% decide not to enter (or more accurately, after entry decision they draw
% their state, and then decide to cancel/abort their entry).
%%

GEPriceParamNames={'p','Ne'}; 

GeneralEqmEqnParamNames(1).Names={'beta'};
GeneralEqmEqn_CondlEntry = @(ValueFn,GEprices,beta) beta*ValueFn-0;
GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.


 heteroagentoptions.specialgeneqmcondn={'condlentry','entry'};
 GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry}; 
  
 %% 
  
%  FnsToEvaluateParamNames(1).Names={'p', 'w','alpha','gamma'};
%FnsToEvaluateFn_1 = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
%    z1_val*(aprime_val^alpha)*...
%    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma); 

%AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
%    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
%    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
%FnsToEvaluate={FnsToEvaluateFn_1};

%GeneralEqmEqnParamNames(3).Names={'p'};
%GeneralEqmEqn_1 = @(AggVars,GEprices,p) 1/AggVars-p;

% heteroagentoptions.specialgeneqmcondn={'condlentry','entry',0};
GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry}; 
%% Find equilibrium prices

n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input 

[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1...
    (n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions,...
    EntryExitParamNames);

Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
Params.ebar=p_eqm.ebar;

%% Value Function, Policy and Firm Distribution in GE
[V,Policy,ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames,ReturnFnParamNames, vfoptions);
Params.oneminuslambda=1-ExitPolicy;
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

% STEP 3 - Now, set K_nfa=K_hh-K_total 
K_nfa=K_hh-K_total;

%% 
%Set r_market back to it's actual value
Params.r_market=Params.r_international;


%%

%
fprintf('\n')
fprintf('Capital Market Outcomes: \n')
fprintf('The interest rate for household is %8.2f, the international interest rate is r%8.2f, and the interest rate on earmarked credit is %8.2f \n', 100*Params.r_hh, 100*Params.r_international, 100*Params.r_ear)
fprintf('Households hold capital of K_hh=%8.4f, while firms use capital of K_total=%8.4f; the difference is made up by net-foreign capital of K_nfa=%8.4f \n',K_hh, K_total, K_nfa)
fprintf('\n')

%%

%figure(2)
%plot(a_grid, cumnansum(nansum(StationaryDist.pdf(:,:,1),2)))
%hold on
%plot(a_grid, cumnansum(nansum(StationaryDist.pdf(:,:,2),2)),'.')
%hold off
%legend('Non-earmarked','Earmarked')
% Also looks a bit stupid with these params and grid, but earmarked are
% ending up with higher capital so that is promising start! :)

%ear =  cumnansum(nansum(StationaryDist.pdf(:,:,2),2));
%nonear =  cumnansum(nansum(StationaryDist.pdf(:,:,1),2));

% Check the average productivity of the earmarked firms 
%AvgProdOfNonEarmarked=nansum(s_grid'.*(nansum(StationaryDist.pdf(:,:,1),1)))
%AvgProdOfEarmarked=nansum(s_grid'.*(nansum(StationaryDist.pdf(:,:,2),1)))

%AvgProdOfNonEarmarked*nonear(end)+ AvgProdOfEarmarked*ear(end)

%figure;
%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(a_grid,cumnansum(nansum(StationaryDist.pdf(:,:,2),2))./ear(end),':','LineWidth',1)
%hold on;
%plot(a_grid,cumnansum(nansum(StationaryDist.pdf(:,:,1),2))./nonear(end),'-','LineWidth',1)
%legend('Earmarked','Non-Earmarked', 'Location', 'Best');



%figure;
%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(s_grid,nansum(StationaryDist.pdf(:,:,2),1),':','LineWidth',1)
%hold on;
%plot(s_grid,nansum(StationaryDist.pdf(:,:,1),1),'-','LineWidth',1)
%xlim([min(s_grid) max(s_grid)])
%legend('Earmarked','Non-Earmarked', 'Location', 'Best');
%% 

% CAPITAL
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass)...
    aprime_val; 

%OUTPUT
FnsToEvaluateParamNames(2).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma);  

%LABOR
FnsToEvaluateParamNames(3).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    ((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));

%% SUB

% CAPITAL WITH SUBSIDY
FnsToEvaluateParamNames(4).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_SUBkbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==1)*aprime_val; 

% OUTPUT WITH SUBSIDY
FnsToEvaluateParamNames(5).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_SUBoutput = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==1)*z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma);

%LABOR WITH SUBSIDY
FnsToEvaluateParamNames(6).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_SUBnbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==1).*(((z1_val*p*gamma)^(1/(1-gamma))) *(aprime_val^(alpha/(1-gamma))));


%% TAX


% CAPITAL WITHOUT SUBSIDY
FnsToEvaluateParamNames(7).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXkbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*aprime_val;

% OUTPUT WITHOUT SUBSIDY
FnsToEvaluateParamNames(8).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXoutput = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*(z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma));

%LABOR WITHOUT SUBSIDY
FnsToEvaluateParamNames(9).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXnbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0).*(((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)));

%% Check

% JUST TO CHECK N FIRMS
FnsToEvaluateParamNames(10).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_num = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*2 + (z2_val==1)*3;

% Subsity costs
FnsToEvaluateParamNames(11).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_cost = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
    (z2_val==1)*(r_market-r_ear)*aprime_val;

% TFP 
FnsToEvaluateParamNames(12).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_tfp = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
    z1_val;

% Subsity costs
FnsToEvaluateParamNames(13).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_eartfp = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
    (z2_val==1)*z1_val;

% TFP 
FnsToEvaluateParamNames(14).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_noneartfp = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
     (z2_val==0)*z1_val;


FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_output, FnsToEvaluateFn_nbar,...
       FnsToEvaluateFn_SUBkbar, FnsToEvaluateFn_SUBoutput, FnsToEvaluateFn_SUBnbar,...
    FnsToEvaluateFn_TAXkbar, FnsToEvaluateFn_TAXoutput, FnsToEvaluateFn_TAXnbar,...
    FnsToEvaluateFn_num,FnsToEvaluateFn_cost, FnsToEvaluateFn_tfp,...
    FnsToEvaluateFn_eartfp,FnsToEvaluateFn_noneartfp};

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(StationaryDist,...
    Policy, FnsToEvaluate, Params,...
    FnsToEvaluateParamNames, n_d, n_a, n_z,...
    [], a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);


ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate,...
    Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,...
    simoptions.parallel,simoptions,EntryExitParamNames);

% kbar output nbar
% Normal 1 2 3
% SUB 4 5 6
% TAX 7 8 9
%% Agggregate Values
%Output.Y=AggVars(2);
%Output.N=AggVars(3);
%Output.K=AggVars(1);
%Output.KdivY=Output.K/Output.Y;
%Output.TFP=(Output.Y/((Output.K^ Params.alpha)*(Output.N^ Params.gamma)));

%% Agggregate Values without Subsidy
%TAX.Output.Y=AggVars(8);
%TAX.Output.N=AggVars(9);
%TAX.Output.K=AggVars(7);
%TAX.Output.KdivY=TAX.Output.K/TAX.Output.Y;
%TAX.Output.TFP=(TAX.Output.Y/((TAX.Output.K^ Params.alpha)*(TAX.Output.N^ Params.gamma)));

%% Agggregate Values with Subsidy
%SUB.Output.Y=AggVars(5);
%SUB.Output.N=AggVars(6);
%SUB.Output.K=AggVars(4);
%SUB.Output.KdivY=SUB.Output.K/SUB.Output.Y;
%SUB.Output.TFP=(SUB.Output.Y/((SUB.Output.K^ Params.alpha)*(SUB.Output.N^ Params.gamma)));

%% Average values

Output.perN=AggVars(3)/StationaryDist.mass;
Output.perK=AggVars(1)/StationaryDist.mass;

firms_sub=100*nansum(nansum(nansum(StationaryDist.pdf(shiftdim(ValuesOnGrid(10,:,:,:),1)==3))));
firms_tax=100*nansum(nansum(nansum(StationaryDist.pdf(shiftdim(ValuesOnGrid(10,:,:,:),1)==2))));

Percentage_tax = [firms_tax  firms_sub  100-firms_tax+firms_sub ] ;

%%

nbarValues=shiftdim(ValuesOnGrid(3,:,:,:),1);
normalize_employment=min(nonzeros(nbarValues)); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=nbarValues./(normalize_employment);

ProbnbarValues=nansum(nansum(shiftdim(ValuesOnGrid(12,:,:,:),1).*...
shiftdim(ProbDensityFns(12,:,:,:),1),3));

% SUB
SUBnbarValues=shiftdim(ValuesOnGrid(6,:,:,:),1);
SUBnbarValues=SUBnbarValues./(normalize_employment);
%SUBnbarValues=SUBnbarValues.*...
%shiftdim(ProbDensityFns(6,:,:,:),1);

% non sub
NONnbarValues=shiftdim(ValuesOnGrid(9,:,:,:),1);
NONnbarValues=NONnbarValues./(normalize_employment);
%NONnbarValues=NONnbarValues.*...
%shiftdim(ProbDensityFns(9,:,:,:),1);
%%
figure;
%subplot(1,2,1)
plot(s_grid,nanmean(nanmean(SUBnbarValues(:,:,:),3)));
xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,nanmean(nanmean(NONnbarValues(:,:,:),3)),'-r');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('earmarked','non-earmarked', 'Location', 'northwest')
%%
Partion1Indicator=logical((nbarValues>0).*(nbarValues<5));
Partion2Indicator=logical((nbarValues>=5).*(nbarValues<50));
Partion3Indicator=logical(nbarValues>=50);


if ((nansum(nansum(nansum(Partion1Indicator+Partion2Indicator+Partion3Indicator)))) - prod(n_z)*(n_a) > 1e-3)
    error('error')
end

ShareOfEstablishments(1)=100*nansum(nansum(nansum(StationaryDist.pdf(Partion1Indicator))))...
    ./nansum(nansum(nansum(StationaryDist.pdf(logical(nbarValues>0)))));
ShareOfEstablishments(2)=100*nansum(nansum(nansum(StationaryDist.pdf(Partion2Indicator))))...
   ./ nansum(nansum(nansum(StationaryDist.pdf(logical(nbarValues>0)))));
ShareOfEstablishments(3)=100*nansum(nansum(nansum(StationaryDist.pdf(Partion3Indicator))))...
    ./nansum(nansum(nansum(StationaryDist.pdf(logical(nbarValues>0)))));
ShareOfEstablishments(4)=100*nansum(nansum(nansum(StationaryDist.pdf(logical(nbarValues>0)))))...
    ./nansum(nansum(nansum(StationaryDist.pdf(logical(nbarValues>0)))));

Output_pdf=shiftdim(ProbDensityFns(2,:,:,:),1);
ShareOfOutput(1)=100*nansum(nansum(nansum(Output_pdf(Partion1Indicator))));
ShareOfOutput(2)=100*nansum(nansum(nansum(Output_pdf(Partion2Indicator))));
ShareOfOutput(3)=100*nansum(nansum(nansum(Output_pdf(Partion3Indicator))));
ShareOfOutput(4)=100*nansum(nansum(nansum(Output_pdf)));

Labour_pdf=shiftdim(ProbDensityFns(3,:,:,:),1);
ShareOfLabour(1)=100*nansum(nansum(nansum(Labour_pdf(Partion1Indicator))));
ShareOfLabour(2)=100*nansum(nansum(nansum(Labour_pdf(Partion2Indicator))));
ShareOfLabour(3)=100*nansum(nansum(nansum(Labour_pdf(Partion3Indicator))));
ShareOfLabour(4)=100*nansum(nansum(nansum(Labour_pdf)));


Capital_pdf=shiftdim(ProbDensityFns(1,:,:,:),1);
ShareOfCapital(1)=100*nansum(nansum(nansum(Capital_pdf(Partion1Indicator))));
ShareOfCapital(2)=100*nansum(nansum(nansum(Capital_pdf(Partion2Indicator))));
ShareOfCapital(3)=100*nansum(nansum(nansum(Capital_pdf(Partion3Indicator))));
ShareOfCapital(4)=100*nansum(nansum(nansum(Capital_pdf)));

AverageEmployment(1)=100*nansum(nansum(nansum(nbarValues(Partion1Indicator).*...
StationaryDist.pdf(Partion1Indicator))))/nansum(nansum(nansum(nbarValues.*...
StationaryDist.pdf)));
AverageEmployment(2)=100*nansum(nansum(nansum(nbarValues(Partion2Indicator).*...
StationaryDist.pdf(Partion2Indicator))))/nansum(nansum(nansum(nbarValues.*...
StationaryDist.pdf)));
AverageEmployment(3)=100*nansum(nansum(nansum(nbarValues(Partion3Indicator).*...
StationaryDist.pdf(Partion3Indicator))))/nansum(nansum(nansum(nbarValues.*...
StationaryDist.pdf)));
AverageEmployment(4)=100*nansum(nansum(nansum(nbarValues.*...
StationaryDist.pdf)))/nansum(nansum(nansum(nbarValues.*...
StationaryDist.pdf)));
%%

%%
TFP_pdf=shiftdim(ValuesOnGrid(12,:,:,:),1);
ShareOfTFP(1)=nansum(TFP_pdf(Partion1Indicator).*(StationaryDist.pdf(Partion1Indicator)/(nansum(nansum(StationaryDist.pdf(Partion1Indicator))))));
ShareOfTFP(2)=nansum(TFP_pdf(Partion2Indicator).*(StationaryDist.pdf(Partion2Indicator)/(nansum(nansum(StationaryDist.pdf(Partion2Indicator))))));
ShareOfTFP(3)=nansum(TFP_pdf(Partion3Indicator).*(StationaryDist.pdf(Partion3Indicator)/(nansum(nansum(StationaryDist.pdf(Partion3Indicator))))));
ShareOfTFP(4)=nansum(nansum(nansum(TFP_pdf(logical(nbarValues>=0)).*StationaryDist.pdf(nbarValues>=0))));


MinOfTFP(1)=nanmin(nanmin(nanmin(nonzeros(TFP_pdf(Partion1Indicator)))));
MinOfTFP(2)=nanmin(nanmin(nanmin(nonzeros(TFP_pdf(Partion2Indicator)))));
MinOfTFP(3)=nanmin(nanmin(nanmin(nonzeros(TFP_pdf(Partion3Indicator)))));



TFP_ear = nansum(nansum(TFP_pdf(:,:,2).*(StationaryDist.pdf(:,:,2)/(nansum(nansum(StationaryDist.pdf(:,:,2)))))));

TFP_nonear = nansum(nansum(TFP_pdf(:,:,1).*(StationaryDist.pdf(:,:,1))/(nansum(nansum(StationaryDist.pdf(:,:,1))))));

%%
SUBStationaryDist.pdf=squeeze(StationaryDist.pdf(:,:,2));
Partion1Indicator=squeeze(Partion1Indicator(:,:,2));
Partion2Indicator=squeeze(Partion2Indicator(:,:,2));
Partion3Indicator=squeeze(Partion3Indicator(:,:,2));

SUBShareOfEstablishments(1)=100*(nansum(nansum(SUBStationaryDist.pdf(Partion1Indicator)))/nansum(nansum(SUBStationaryDist.pdf)));
SUBShareOfEstablishments(2)=100*(nansum(nansum(SUBStationaryDist.pdf(Partion2Indicator)))/nansum(nansum(SUBStationaryDist.pdf)));
SUBShareOfEstablishments(3)=100*(nansum(nansum(SUBStationaryDist.pdf(Partion3Indicator)))/nansum(nansum(SUBStationaryDist.pdf)));
SUBShareOfEstablishments(4)=100*(nansum(nansum(SUBStationaryDist.pdf))/nansum(nansum(SUBStationaryDist.pdf)));

SUBOutput_pdf=shiftdim(ProbDensityFns(2,:,:,2),1);
SUBShareOfOutput(1)=100*(nansum(nansum(SUBOutput_pdf(Partion1Indicator)))/nansum(nansum(SUBOutput_pdf)));
SUBShareOfOutput(2)=100*(nansum(nansum(SUBOutput_pdf(Partion2Indicator)))/nansum(nansum(SUBOutput_pdf)));
SUBShareOfOutput(3)=100*(nansum(nansum(SUBOutput_pdf(Partion3Indicator)))/nansum(nansum(SUBOutput_pdf)));
SUBShareOfOutput(4)=100*(nansum(nansum(SUBOutput_pdf))/nansum(nansum(SUBOutput_pdf)));

SUBLabour_pdf=shiftdim(ProbDensityFns(3,:,:,2),1);
SUBShareOfLabour(1)=100*(nansum(nansum(SUBLabour_pdf(Partion1Indicator)))/nansum(nansum(SUBLabour_pdf)));
SUBShareOfLabour(2)=100*(nansum(nansum(SUBLabour_pdf(Partion2Indicator)))/nansum(nansum(SUBLabour_pdf)));
SUBShareOfLabour(3)=100*(nansum(nansum(SUBLabour_pdf(Partion3Indicator)))/nansum(nansum(SUBLabour_pdf)));
SUBShareOfLabour(4)=100*(nansum(nansum(SUBLabour_pdf))/nansum(nansum(SUBLabour_pdf)));


SUBCapital_pdf=shiftdim(ProbDensityFns(1,:,:,2),1);
SUBShareOfCapital(1)=100*(nansum(nansum(SUBCapital_pdf(Partion1Indicator)))/(nansum(nansum(SUBCapital_pdf))));
SUBShareOfCapital(2)=100*(nansum(nansum(SUBCapital_pdf(Partion2Indicator)))/(nansum(nansum(SUBCapital_pdf))));
SUBShareOfCapital(3)=100*(nansum(nansum(SUBCapital_pdf(Partion3Indicator)))/(nansum(nansum(SUBCapital_pdf))));


%%
%fprintf('Distribution statistics of benchmark economy  \n');
%fprintf('                               <5     5 to 49     >=50   total\n');
%fprintf('Share of establishments  %8.2f  %8.2f  %8.2f  %8.2f  \n', ShareOfEstablishments);
%fprintf('Share of output          %8.2f  %8.2f  %8.2f  %8.2f\n', ShareOfOutput);
%fprintf('Share of labor          %8.2f  %8.2f  %8.2f  %8.2f\n', ShareOfLabour);
%fprintf('Share of capital         %8.2f  %8.2f  %8.2f  %8.2f\n', ShareOfCapital);
%fprintf('Share of employment      %8.2f  %8.2f  %8.2f  %8.2f\n', AverageEmployment);

%% Display some output about the solution

%fprintf('The equilibrium output price is p=%.4f \n', Params.p)
%fprintf('The equilibrium value for the mass of entrants is Ne=%.4f \n', Params.Ne)

%fprintf('Total Output is Y=%.4f \n', Output.Y)
%fprintf('Labor is n=%.4f \n', Output.N)
%fprintf('Capital is k=%.4f \n', Output.K)
%fprintf('Total Factor Productivity is TFP=%.4f \n', Output.TFP)


%fprintf('Percentage of firms with\n')
%fprintf('Market Rate     Subsidized Rate\n')
%fprintf('%9.2f  %12.2f   \n',Percentage_tax )

%fprintf('   Total(just for checking)  \n' )
%fprintf('%9.2f    \n', nansum(Percentage_tax))


%fprintf('                     Total  with r_ear   with r_market\n');
%fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[Output.TFP SUB.Output.TFP TAX.Output.TFP]);
%fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [Output.Y SUB.Output.Y TAX.Output.Y]);
%fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [Output.N SUB.Output.N TAX.Output.N]);
%fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [Output.K SUB.Output.K TAX.Output.K]);



%hold on;
%plot(s_grid,(nansum(squeeze(StationaryDist.pdf(:,:,1)),1)),'b')
%hold on;
%plot(s_grid,(nansum(squeeze(StationaryDist.pdf(:,:,2)),1)),':')

lambda = nansum(nansum(nansum(Params.oneminuslambda(:,:,:)==1)))/nansum(nansum(nansum(Params.oneminuslambda(:,:,:)>=0)));
probenter = nansum(nansum(nansum(Params.ebar(:,:,:)==1)))/nansum(nansum(nansum(Params.ebar(:,:,:)>=0)));


%toc;
