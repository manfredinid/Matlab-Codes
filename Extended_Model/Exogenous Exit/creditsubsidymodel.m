
%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% Model presented in he TBrazilian Slump and the Government-driven 
% Credit Expansion (2020)

%% Initial setups
%clear all;
%close all;
%Parallel=1; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
%tic;

%% Toolkit options 

vfoptions.parallel=Parallel;
simoptions.parallel=Parallel;
heteroagentoptions.verbose=1;
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;

%% Parameters

% Preferences 
Params.beta=0.9798;% Discount rate

% Firm-level technology
Params.alpha=0.399;  % Capital share
Params.gamma=0.491; % alpha + gamma must be ~= 1
Params.delta=0.025; % Depreciation rate of physical capital
Params.cf=0.2;%0.0392; % Fixed cost of production

% Adjustment cost of capital
Params.adjustcostparam = 3.219;

% Entry and Exit
Params.ce=1.05;%.00991; % Fixed cost of entry 
Params.lambda= 1-(1-0.1859)^(1/4); % Probability of firm exit
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

n_s=13;
n_a=501;
% n_psi is two since psi \in {0,1}

%% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

%Params.r_ear=0.02; % Interest rate on earmarked credit
%Params.g_ear=0.4; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.

%% Productivity (s)

rhoeps = 0.9; % persistence
evallowpareto = 1.5; % lower bound
evalhighpareto = 2.75;%upper bound
eparampareto = 5.8;% shape parameter
% lower eparampreto -- less small firms
s_grid = linspace(evallowpareto,evalhighpareto,n_s);
rand('state',1)
[pistar_s, pi_s]= paretojo(n_s, s_grid, eparampareto, rhoeps);
pistar_s=pistar_s';
s_grid=s_grid';

% Earmarked credit grid
psi_grid=[0;1]; % Using this as a {0,1} helps, e.g., add up earmarked capital.
pi_psi=[1,0;0,1];

%figure; plot(s_grid, pistar_s)

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
k_high = log(((Params.alpha/(Params.r_ear))^((1-Params.gamma)/(1-Params.gamma-Params.alpha)))*(Params.gamma^((Params.gamma)/(1-Params.gamma-Params.alpha)))*(s_grid(end)^(1/(1-Params.gamma-Params.alpha))))/log(10);
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

%% normalization
Params.w=1; % Normalization


%% Aspects of the Endogenous entry
% Exit is exogenous with probability lambda
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;

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

%%


%% Return Function
ReturnFn=@(kprime_val, k_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma,delta, cf, adjustcostparam) ExistingFirm_ReturnFn(kprime_val, k_val,s_val, psi_val, p,w,r_market,r_ear, alpha,gamma,delta, cf, adjustcostparam);
ReturnFnParamNames={'p','w','r_market','r_ear', 'alpha','gamma','delta', 'cf', 'adjustcostparam'}; %It is important that these are in same order as they appear in 'ExistingFirm_ReturnFn'

%% General Equilibrium Equations
%Now define the functions for the General Equilibrium conditions

GEPriceParamNames={'ebar'};

FnsToEvaluateParamNames(1).Names={'p','alpha','gamma'};
FnsToEvaluateFn_nbar1 = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,alpha,gamma)...
((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)); 
FnsToEvaluate={FnsToEvaluateFn_nbar1};

GEPriceParamNames={'Ne','p'};

% Conditional entry
GeneralEqmEqnParamNames(1).Names={'beta'};
GeneralEqmEqn_CondlEntry = @(ValueFn,GEprices,beta) beta*ValueFn-0;

% Entry
GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,beta,ce) beta*EValueFn-ce; 

% Labor Market Equilibrium
GeneralEqmEqnParamNames(3).Names={};
GeneralEqmEqn_LabourMarket = @(AggVars,GEprices) 1-AggVars(1);

heteroagentoptions.specialgeneqmcondn={'condlentry','entry',0};

GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry,GeneralEqmEqn_LabourMarket};


%% Find equilibrium prices

heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input 

[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1...
    (0, n_a, n_z, n_p, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params,...
    DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames,...
    GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions,...
    simoptions, vfoptions, EntryExitParamNames);
%%
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
Params.ebar=p_eqm.ebar;

%[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1...
%    (0, n_a, n_z, n_p, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params,...
%    DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames,...
%    GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions,...
%    simoptions, vfoptions, EntryExitParamNames);
%%
%Params.p=p_eqm.p;
%Params.Ne=gather(p_eqm.Ne);
%Params.ebar=p_eqm.ebar;

%% Value Function, Policy and Firm Distribution in GE
[V,Policy]=ValueFnIter_Case1(0,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn,...
    Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
StationaryDist=StationaryDist_Case1(Policy,0,n_a,n_z,pi_z,...
    simoptions, Params, EntryExitParamNames);

%%
FnsToEvaluateParamNames(1).Names={'p','alpha','gamma'};
FnsToEvaluateFn_nbar1 = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,alpha,gamma)...
((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)); 
FnsToEvaluate={FnsToEvaluateFn_nbar1};
%%
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid,z_grid, simoptions.parallel,simoptions,EntryExitParamNames)



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
%plot(a_grid, cumsum(sum(StationaryDist.pdf(:,:,1),2)))
%hold on
%plot(a_grid, cumsum(sum(StationaryDist.pdf(:,:,2),2)),'.')
%hold off
%legend('Non-earmarked','Earmarked')
% Also looks a bit stupid with these params and grid, but earmarked are
% ending up with higher capital so that is promising start! :)

%ear =  cumsum(sum(StationaryDist.pdf(:,:,2),2));
%nonear =  cumsum(sum(StationaryDist.pdf(:,:,1),2));

% Check the average productivity of the earmarked firms 
%AvgProdOfNonEarmarked=sum(s_grid'.*(sum(StationaryDist.pdf(:,:,1),1)))
%AvgProdOfEarmarked=sum(s_grid'.*(sum(StationaryDist.pdf(:,:,2),1)))

%AvgProdOfNonEarmarked*nonear(end)+ AvgProdOfEarmarked*ear(end)

%figure;
%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(a_grid,cumsum(sum(StationaryDist.pdf(:,:,2),2))./ear(end),':','LineWidth',1)
%hold on;
%plot(a_grid,cumsum(sum(StationaryDist.pdf(:,:,1),2))./nonear(end),'-','LineWidth',1)
%legend('Earmarked','Non-Earmarked', 'Location', 'Best');



%figure;
%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(s_grid,sum(StationaryDist.pdf(:,:,2),1),':','LineWidth',1)
%hold on;
%plot(s_grid,sum(StationaryDist.pdf(:,:,1),1),'-','LineWidth',1)
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
    (z2_val==1)*((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));


%% TAX


% CAPITAL WITHOUT SUBSIDY
FnsToEvaluateParamNames(7).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXkbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*aprime_val;

% OUTPUT WITHOUT SUBSIDY
FnsToEvaluateParamNames(8).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXoutput = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma);

%LABOR WITHOUT SUBSIDY
FnsToEvaluateParamNames(9).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXnbar = @(aprime_val,a_val,z1_val,z2_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));

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
Output.Y=AggVars(2);
Output.N=AggVars(3);
Output.K=AggVars(1);
Output.KdivY=Output.K/Output.Y;
Output.TFP=(Output.Y/((Output.K^ Params.alpha)*(Output.N^ Params.gamma)));

%% Agggregate Values without Subsidy
TAX.Output.Y=AggVars(8);
TAX.Output.N=AggVars(9);
TAX.Output.K=AggVars(7);
TAX.Output.KdivY=TAX.Output.K/TAX.Output.Y;
TAX.Output.TFP=(TAX.Output.Y/((TAX.Output.K^ Params.alpha)*(TAX.Output.N^ Params.gamma)));

%% Agggregate Values with Subsidy
SUB.Output.Y=AggVars(5);
SUB.Output.N=AggVars(6);
SUB.Output.K=AggVars(4);
SUB.Output.KdivY=SUB.Output.K/SUB.Output.Y;
SUB.Output.TFP=(SUB.Output.Y/((SUB.Output.K^ Params.alpha)*(SUB.Output.N^ Params.gamma)));

%% Average values

Output.perN=AggVars(3)/StationaryDist.mass;
Output.perK=AggVars(1)/StationaryDist.mass;

firms_sub=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(10,:,:,:),1)==3))));
firms_tax=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(10,:,:,:),1)==2))));

Percentage_tax = [firms_tax  firms_sub] ;

%%

nbarValues=shiftdim(ValuesOnGrid(3,:,:,:),1);
normalize_employment=nanmin(nonzeros(nbarValues)); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
if normalize_employment<1
nbarValues=nbarValues./(normalize_employment);
end


ProbnbarValues=sum(sum(shiftdim(ValuesOnGrid(12,:,:,:),1).*...
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
plot(s_grid,nanmean(nanmean(SUBnbarValues(:,:,:),3)),':k');
%xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,nanmean(nanmean(NONnbarValues(:,:,:),3)),'-k');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('earmarked','non-earmarked', 'Location', 'northwest')
%%




Partion1Indicator=logical(nbarValues<5);
Partion2Indicator=logical((nbarValues>=5).*(nbarValues<50));
Partion3Indicator=logical(nbarValues>=50);


if ((sum(sum(sum(Partion1Indicator+Partion2Indicator+Partion3Indicator)))) - prod(n_z)*(n_a) > 1e-3)
    error('error')
end

ShareOfEstablishments(1)=100*sum(sum(sum(StationaryDist.pdf(Partion1Indicator))));
ShareOfEstablishments(2)=100*sum(sum(sum(StationaryDist.pdf(Partion2Indicator))));
ShareOfEstablishments(3)=100*sum(sum(sum(StationaryDist.pdf(Partion3Indicator))));
ShareOfEstablishments(4)=100*sum(sum(sum(StationaryDist.pdf)));

Output_pdf=shiftdim(ProbDensityFns(2,:,:,:),1);
ShareOfOutput(1)=100*sum(sum(sum(Output_pdf(Partion1Indicator))));
ShareOfOutput(2)=100*sum(sum(sum(Output_pdf(Partion2Indicator))));
ShareOfOutput(3)=100*sum(sum(sum(Output_pdf(Partion3Indicator))));
ShareOfOutput(4)=100*sum(sum(sum(Output_pdf)));

Labour_pdf=shiftdim(ProbDensityFns(3,:,:,:),1);
ShareOfLabour(1)=100*sum(sum(sum(Labour_pdf(Partion1Indicator))));
ShareOfLabour(2)=100*sum(sum(sum(Labour_pdf(Partion2Indicator))));
ShareOfLabour(3)=100*sum(sum(sum(Labour_pdf(Partion3Indicator))));
ShareOfLabour(4)=100*sum(sum(sum(Labour_pdf)));



Capital_pdf=shiftdim(ProbDensityFns(1,:,:,:),1);
ShareOfCapital(1)=100*sum(sum(sum(Capital_pdf(Partion1Indicator))));
ShareOfCapital(2)=100*sum(sum(sum(Capital_pdf(Partion2Indicator))));
ShareOfCapital(3)=100*sum(sum(sum(Capital_pdf(Partion3Indicator))));
ShareOfCapital(4)=100*sum(sum(sum(Capital_pdf)));

AverageEmployment(1)=100*sum(sum(sum(nbarValues(Partion1Indicator).*...
StationaryDist.pdf(Partion1Indicator))))/sum(sum(sum(nbarValues.*...
StationaryDist.pdf)));
AverageEmployment(2)=100*sum(sum(sum(nbarValues(Partion2Indicator).*...
StationaryDist.pdf(Partion2Indicator))))/sum(sum(sum(nbarValues.*...
StationaryDist.pdf)));
AverageEmployment(3)=100*sum(sum(sum(nbarValues(Partion3Indicator).*...
StationaryDist.pdf(Partion3Indicator))))/sum(sum(sum(nbarValues.*...
StationaryDist.pdf)));
AverageEmployment(4)=100*sum(sum(sum(nbarValues.*...
StationaryDist.pdf)))/sum(sum(sum(nbarValues.*...
StationaryDist.pdf)));


TFP_pdf=shiftdim(ValuesOnGrid(12,:,:,:),1);
ShareOfTFP(1)=sum(TFP_pdf(Partion1Indicator).*(StationaryDist.pdf(Partion1Indicator)/(sum(sum(StationaryDist.pdf(Partion1Indicator))))));
ShareOfTFP(2)=sum(TFP_pdf(Partion2Indicator).*(StationaryDist.pdf(Partion2Indicator)/(sum(sum(StationaryDist.pdf(Partion2Indicator))))));
ShareOfTFP(3)=sum(TFP_pdf(Partion3Indicator).*(StationaryDist.pdf(Partion3Indicator)/(sum(sum(StationaryDist.pdf(Partion3Indicator))))));
%ShareOfTFP(4)=sum(sum(sum(TFP_pdf.*StationaryDist.pdf)));
ShareOfTFP(4)=sum(sum(sum(TFP_pdf.*(shiftdim(ProbDensityFns(12,:,:,:),1)))));


TFP_ear = sum(sum(TFP_pdf(:,:,2).*(StationaryDist.pdf(:,:,2)/(sum(sum(StationaryDist.pdf(:,:,2)))))));

TFP_nonear = sum(sum(TFP_pdf(:,:,1).*(StationaryDist.pdf(:,:,1)/(sum(sum(StationaryDist.pdf(:,:,1)))))));

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
%fprintf('%9.2f    \n', sum(Percentage_tax))


%fprintf('                     Total  with r_ear   with r_market\n');
%fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[Output.TFP SUB.Output.TFP TAX.Output.TFP]);
%fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [Output.Y SUB.Output.Y TAX.Output.Y]);
%fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [Output.N SUB.Output.N TAX.Output.N]);
%fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [Output.K SUB.Output.K TAX.Output.K]);



%hold on;
%plot(s_grid,(sum(squeeze(StationaryDist.pdf(:,:,1)),1)),'b')
%hold on;
%plot(s_grid,(sum(squeeze(StationaryDist.pdf(:,:,2)),1)),':')


%toc;
