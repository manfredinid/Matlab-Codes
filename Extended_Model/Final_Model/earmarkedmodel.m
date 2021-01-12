%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% The model has 3 exogenous states and three types of exit

%% Toolkit options 

vfoptions.agententryandexit=1;
vfoptions.endogenousexit=2;
vfoptions.parallel=Parallel;
simoptions.parallel=Parallel;
simoptions.agententryandexit=1;
simoptions.endogenousexit=2;
simoptions.iterate=1;

%% Parameters

% Preferences 
Params.beta=0.9798;% Discount rate

% Firm-level technology
Params.alpha=0.399;  % Capital share
Params.gamma=0.491; % alpha + gama must be ~= 1
Params.delta=0.025; % Depreciation rate of physical capital
Params.cf=0.05; % Fixed cost of production

Params.w=1; % Normalization

% Adjustment cost of capital
Params.adjustcostparam = 3.219;

% Entry and Exit
Params.ce=2; % Fixed cost of entry 
Params.phi=0.02; % Continuation fixed cost for firms facing endogenous exit decision




%% States

% The model has three states, one endogenous state (capital), and two
% exogenous states (productivity and subsidies)

n_s=10;
n_a=280;
% n_psi is two since psi \in {0,1}
% n_tau is two since psi \in {0,1}

%% Grids and Transition Matrices

% Productivity (s)
% Exogenous AR(1) process on (log) productivity
% logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)

rhoeps = 0.95; % persistence
evallowpareto = 1.8; % lower bound
evalhighpareto = 2.75;%upper bound
eparampareto = 5.9;% shape parameter
% lower eparampreto -- less small firms
s_grid = linspace(evallowpareto,evalhighpareto,n_s);
rand('state',1)
[pistar_s, pi_s]= paretodistribution(n_s, s_grid, eparampareto, rhoeps);
pistar_s=pistar_s';
s_grid=s_grid';

% Earmarked credit grid
psi_grid=[0;1]; %Using this as a {0,1} helps,e.g., add up earmarked capital
pi_psi=[1,0;0,1];

tau_grid=[0;1]; %Using this as a {0,1} helps,e.g., add up poor credit access
pi_tau=[1,0;0,1];

% Exogenous states (matrix z)
% Transition matrix 
n_z=[n_s,length(psi_grid), length(tau_grid)];
z_grid=[s_grid; psi_grid; tau_grid];

% Transition matrix for the exogenous states
pi_z=kron(pi_tau,kron(pi_psi, pi_s));

% grid for capital
a_grid = [0 logspace(0.0001,6.28,n_a-1)]'; 

% Decision variables
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
%Params.r_international= Main file

%% 3.Market interest rate
Params.r_market=Params.r_international;

%% Potential New Entrants Distribution over the states (s, psi)


pistar_s=ones(size(s_grid))/n_s; % Initial guess
dist=1;
while dist>10^(-7)
    pistar_s_old=pistar_s;
    pistar_s=(pi_s)'*pistar_s;
    dist=max(abs(pistar_s-pistar_s_old));
end

if Parallel==2
    Params.upsilon = zeros([n_a, n_z],'gpuArray');
else
    Params.upsilon = zeros([n_a, n_z]);
end

% Transition matrix 
n_z=[n_s,length(psi_grid), length(tau_grid)];
z_grid=[s_grid; psi_grid; tau_grid];

% Transition matrix for the exogenous states
pi_z=kron(pi_tau,kron(pi_psi, pi_s));

% grid for capital
a_grid = [0 logspace(0.0001,6.28,n_a-1)]'; 

% Decision variables
%There is no d variable
d_grid=[]; 
n_d=0;

A=[Params.g_tau 1-Params.g_tau; 1-Params.g_tau Params.g_tau].*[1-...
    Params.g_ear, Params.g_ear];

B = kron(A,pistar_s);

Params.upsilon(1,:,:,:) = reshape(B,[10,2,2]);


%%% Revisar se eu tenho tipi 1, 2 e 3 so mudando gtau
% FIX FIX FIX

%% Aspects of Entry and Exit 
% Discount Factor
DiscountFactorParamNames={'beta'};

% Exit status
Params.lambda_phi=0.04;    %endogenous exit decision
Params.lambda_infty=0.02; %exogenous exit decision

vfoptions.exitprobabilities={'lambda_phi','lambda_infty'};
simoptions.exitprobabilities=vfoptions.exitprobabilities;

% For firms facing 'endogenous exit decision', there is a continuation cost:

vfoptions.endogenousexitcontinuationcost={'phi'};

% Return Function for Remaning Firms
ReturnFn=@(kprime_val, k_val,s_val, psi_val,tau_val,  p,w,r_market,r_ear,...
    alpha,gamma,delta, ctau, adjustcostparam, cf)RemaningFirm_ReturnFn(kprime_val, k_val,s_val, psi_val, tau_val,  p,w,r_market,r_ear, alpha,gamma, delta, ctau, adjustcostparam, cf);
ReturnFnParamNames={'p','w','r_market','r_ear', 'alpha','gamma','delta','ctau',...
     'adjustcostparam', 'cf'}; 
%It is important that these are in same order as they appear in 'ExistingFirm_ReturnFn'


% Return Function for Exiting Firms
vfoptions.ReturnToExitFn=@(kprime_val, k_val,s_val, psi_val,tau_val,  p,w,r_market,r_ear,...
    alpha,gamma,delta, ctau, adjustcostparam, cf) ExistingFirm_ReturnFn(kprime_val, k_val,s_val, psi_val,tau_val,  p,w,r_market,r_ear,...
    alpha,gamma,delta, ctau, adjustcostparam, cf);
vfoptions.ReturnToExitFnParamNames={'p','w','r_market','r_ear', 'alpha','gamma','delta','ctau',...
     'adjustcostparam', 'cf'};

%% Check

[V, Policy, PolicyWhenExiting, ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.notexit=1-ExitPolicy;



%% Set up variables for entry and exit

% Probability of being in the (k, s, psi, tau) category
EntryExitParamNames.DistOfNewAgents={'upsilon'};

% Exit decision
EntryExitParamNames.CondlProbOfSurvival={'notexit'};
Params.notexit=1-ExitPolicy;

% Conditional entry will allow for different productivity cutoffs for 
%new entrants depending on earmarked-vs-nonearmarked.
EntryExitParamNames.CondlEntryDecisions={'ebar'};
% Takes value of one for enter, zero for not-enter. This is just an initial
%guess as the actual decisions are determined as part of general equilibrium.
Params.ebar=ones([n_a,n_z],'gpuArray'); 

% Mass of new entrants
EntryExitParamNames.MassOfNewAgents={'Ne'};


% Some checks
disp('upsilon size')
disp(size(Params.upsilon))
disp('sum of upsilon')
disp(sum(Params.upsilon(:)))


%% General Equilibrium Equations
%Now define the functions for the General Equilibrium conditions

GEPriceParamNames={'ebar'};
FnsToEvaluateParamNames(1).ExitStatus=[1,1,1,1];
FnsToEvaluateParamNames(1).Names={'p','alpha','gamma'};
FnsToEvaluateFn_nbar1 = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,alpha,gamma)...
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
GeneralEqmEqn_LabourMarket = @(AggVars,GEprices) 1-AggVars;

heteroagentoptions.specialgeneqmcondn={'condlentry','entry',0};

GeneralEqmEqns={GeneralEqmEqn_CondlEntry,GeneralEqmEqn_Entry,GeneralEqmEqn_LabourMarket};

%% Check

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);
sum(sum(sum(StationaryDist.pdf(:,:,:,1))))
sum(sum(sum(StationaryDist.pdf(:,:,1,:))))

sum(sum(sum(StationaryDist.pdf(:,:,1,1))))
%% Find equilibrium prices
heteroagentoptions.verbose=1;
simoptions.endogenousexit=2;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% NOTE: EntryExitParamNames has to be passed as an additional input 
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);

Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
Params.ebar=p_eqm.ebar;

%% Value Function, Policy and Firm Distribution in GE
[V, Policy, PolicyWhenExiting, ExitPolicy]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
Params.notexit=1-ExitPolicy;

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions,Params,EntryExitParamNames);

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

% Find Aggregate Capital with r = r_market
%clear FnsToEvaluateParamNames
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_Kfirm = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass)...
    aprime_val; 
FnsToEvaluate={FnsToEvaluateFn_Kfirm};


K_total=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, Parallel, simoptions, EntryExitParamNames, PolicyWhenExiting);

%%        
% STEP 2 - Let K_HH=K_firm(r_HH)  
% Params.rhhminusdelta=1/Params.beta-1, that is general eqm result in 
% complete market models
% Params.r_hh=Params.rhhminusdelta+Params.delta;
Params.r_market=Params.r_hh;
% Policy function with r = r_HH
[~,Policy_rhh,~,~]=ValueFnIter_Case1(n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);
K_hh=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy_rhh, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, Parallel, simoptions, EntryExitParamNames, PolicyWhenExiting);

% STEP 3 - Now, set K_nfa=K_hh-K_total 
K_nfa=K_hh-K_total;

%Set r_market back to it's actual value
Params.r_market=Params.r_international;

%% Calculating the Model Results

% CAPITAL
FnsToEvaluateParamNames(1).ExitStatus=[1,1,1,1]; % [NoExit, EndogExit-choosenottoexit, EndogExit-choosetoexit,ExoExit] (1 means include 'this' status in calculation, 0 not include) (the default when not specified is [1,1,1,1])
FnsToEvaluateParamNames(1).Names={};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass)...
    aprime_val; 

%OUTPUT
FnsToEvaluateParamNames(2).ExitStatus=[1,1,1,1];
FnsToEvaluateParamNames(2).Names={'p','alpha','gamma'};
FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,alpha,gamma)...
    z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma);  

%LABOR
FnsToEvaluateParamNames(3).ExitStatus=[1,1,1,1];
FnsToEvaluateParamNames(3).Names={'p','alpha','gamma'};
FnsToEvaluateFn_nbar =  @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,alpha,gamma)...
((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)); 

%% SUB

% CAPITAL WITH SUBSIDY
FnsToEvaluateParamNames(4).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_SUBkbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==1)*aprime_val; 

% OUTPUT WITH SUBSIDY
FnsToEvaluateParamNames(5).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_SUBoutput = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==1)*z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma);

%LABOR WITH SUBSIDY
FnsToEvaluateParamNames(6).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_SUBnbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==1).*(((z1_val*p*gamma)^(1/(1-gamma))) *(aprime_val^(alpha/(1-gamma))));


%% TAX


% CAPITAL WITHOUT SUBSIDY
FnsToEvaluateParamNames(7).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXkbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*aprime_val;

% OUTPUT WITHOUT SUBSIDY
FnsToEvaluateParamNames(8).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXoutput = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*(z1_val*(aprime_val^alpha)*...
    (((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma))^gamma));

%LABOR WITHOUT SUBSIDY
FnsToEvaluateParamNames(9).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_TAXnbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0).*(((z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)));

%% Check

% JUST TO CHECK N FIRMS
FnsToEvaluateParamNames(10).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_num = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z2_val==0)*2 + (z2_val==1)*3;

% Subsity costs
FnsToEvaluateParamNames(11).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_cost = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
    (z2_val==1)*(r_market-r_ear)*aprime_val;

% TFP 
FnsToEvaluateParamNames(12).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_tfp = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
    z1_val;

% Subsity costs
FnsToEvaluateParamNames(13).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_eartfp = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
    (z2_val==1)*z1_val;

% TFP 
FnsToEvaluateParamNames(14).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_noneartfp = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
     (z2_val==0)*z1_val;
 
 % TFP 
FnsToEvaluateParamNames(15).Names={'p', 'w','r_market','r_ear','alpha','gamma'};
FnsToEvaluateFn_poor = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,r_market,r_ear,alpha,gamma)...
     (z3_val==1);
 
 %LABOR WITH SUBSIDY
FnsToEvaluateParamNames(16).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_POORnbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z3_val==1).*(((z1_val*p*gamma)^(1/(1-gamma))) *(aprime_val^(alpha/(1-gamma))));


%LABOR WITH SUBSIDY
FnsToEvaluateParamNames(17).Names={'p', 'w','alpha','gamma'};
FnsToEvaluateFn_GOODnbar = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,alpha,gamma)...
    (z3_val==0).*(((z1_val*p*gamma)^(1/(1-gamma))) *(aprime_val^(alpha/(1-gamma))));


FnsToEvaluateParamNames(18).Names={'p', 'w','r_market','r_ear', 'ctau'};
FnsToEvaluateFn_r = @(aprime_val,a_val,z1_val,z2_val,z3_val,AgentDistMass,p,w,r_market,r_ear, ctau)...
    z3_val*ctau + r_ear*z2_val + r_market*(1-z2_val);



FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_output, FnsToEvaluateFn_nbar,...
       FnsToEvaluateFn_SUBkbar, FnsToEvaluateFn_SUBoutput, FnsToEvaluateFn_SUBnbar,...
    FnsToEvaluateFn_TAXkbar, FnsToEvaluateFn_TAXoutput, FnsToEvaluateFn_TAXnbar,...
    FnsToEvaluateFn_num,FnsToEvaluateFn_cost, FnsToEvaluateFn_tfp,...
    FnsToEvaluateFn_eartfp,FnsToEvaluateFn_noneartfp,FnsToEvaluateFn_poor,...
    FnsToEvaluateFn_POORnbar, FnsToEvaluateFn_GOODnbar,FnsToEvaluateFn_r};

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, Parallel, simoptions, EntryExitParamNames, PolicyWhenExiting);


ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1(StationaryDist,...
    Policy, FnsToEvaluate, Params,...
    FnsToEvaluateParamNames, n_d, n_a, n_z,...
    [], a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);


ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate,...
    Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,...
    simoptions.parallel,simoptions,EntryExitParamNames);



%% Average values

Output.perN=AggVars(3)/StationaryDist.mass;
Output.perK=AggVars(1)/StationaryDist.mass;

firms_sub=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(10,:,:,:),1)==3))));
firms_tax=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(10,:,:,:),1)==2))));

Percentage_tax = [firms_tax  firms_sub  firms_tax+firms_sub ] ;

MassOfExitingFirms=sum(sum(sum(StationaryDist.pdf(logical(ExitPolicy)))))*StationaryDist.mass;
ExitRateOfFirms=(MassOfExitingFirms/StationaryDist.mass).*Params.lambda_phi+Params.lambda_infty;


%%

nbarValues=shiftdim(ValuesOnGrid(3,:,:,:,:),1);
normalize_employment=nanmin(nonzeros(nbarValues)); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=nbarValues./(normalize_employment);

% SUB
SUBnbarValues=shiftdim(ValuesOnGrid(6,:,:,:,:),1);
SUBnbarValues=SUBnbarValues./(normalize_employment);

% non sub
NONnbarValues=shiftdim(ValuesOnGrid(9,:,:,:,:),1);
NONnbarValues=NONnbarValues./(normalize_employment);


%%
Partion1Indicator=logical(nbarValues<5);
Partion2Indicator=logical((nbarValues>=5).*(nbarValues<50));
Partion3Indicator=logical(nbarValues>=50);


ShareOfEstablishments(1)=100*sum(sum(sum(sum(StationaryDist.pdf(logical(nbarValues<5))))));
ShareOfEstablishments(2)=100*sum(sum(sum(sum(StationaryDist.pdf(Partion2Indicator)))));
ShareOfEstablishments(3)=100*sum(sum(sum(sum(StationaryDist.pdf(Partion3Indicator)))));
ShareOfEstablishments(4)=ShareOfEstablishments(1)+ShareOfEstablishments(2)+ShareOfEstablishments(3);

Output_pdf=shiftdim(ProbDensityFns(2,:,:,:),1);
ShareOfOutput(1)=100*sum(sum(sum(sum(Output_pdf(logical((nbarValues>0).*(nbarValues<5)))))));
ShareOfOutput(2)=100*sum(sum(sum(sum(Output_pdf(logical((nbarValues>=5).*(nbarValues<50)))))));
ShareOfOutput(3)=100*sum(sum(sum(sum(Output_pdf(logical(nbarValues>=50))))));
ShareOfOutput(4)=ShareOfOutput(1)+ShareOfOutput(2)+ShareOfOutput(3);

Labour_pdf=shiftdim(ProbDensityFns(3,:,:,:),1);
ShareOfLabour(1)=100*sum(sum(sum(sum(Labour_pdf(logical((nbarValues>0).*(nbarValues<5)))))));
ShareOfLabour(2)=100*sum(sum(sum(sum(Labour_pdf(logical((nbarValues>=5).*(nbarValues<50)))))));
ShareOfLabour(3)=100*sum(sum(sum(sum(Labour_pdf(logical(nbarValues>=50))))));
ShareOfLabour(4)=ShareOfLabour(1)+ShareOfLabour(2)+ShareOfLabour(3);


Capital_pdf=shiftdim(ProbDensityFns(1,:,:,:),1);
ShareOfCapital(1)=100*sum(sum(sum(sum(Capital_pdf(logical((nbarValues>0).*(nbarValues<5)))))));
ShareOfCapital(2)=100*sum(sum(sum(sum(Capital_pdf(logical((nbarValues>=5).*(nbarValues<50)))))));
ShareOfCapital(3)=100*sum(sum(sum(sum(Capital_pdf(logical(nbarValues>=50))))));
ShareOfCapital(4)=ShareOfCapital(1)+ShareOfCapital(2)+ShareOfCapital(3);

AverageEmployment(1)=nansum(nansum(nansum(nansum(nbarValues(logical((nbarValues>0).*(nbarValues<5))).*StationaryDist.pdf(logical((nbarValues>0).*(nbarValues<5)))))))/sum(sum(sum(sum(StationaryDist.pdf(logical((nbarValues>0).*(nbarValues<5)))))));
AverageEmployment(2)=nansum(nansum(nansum(nansum(nbarValues(logical((nbarValues>=5).*(nbarValues<50))).*StationaryDist.pdf(logical((nbarValues>=5).*(nbarValues<50)))))))/sum(sum(sum(sum(StationaryDist.pdf(logical((nbarValues>=5).*(nbarValues<50)))))));
AverageEmployment(3)=nansum(nansum(nansum(nansum(nbarValues(logical(nbarValues>=50)).*StationaryDist.pdf(logical(nbarValues>=50))))))/sum(sum(sum(sum(StationaryDist.pdf(logical(nbarValues>=50))))));
AverageEmployment(4)=nansum(nansum(nansum(nansum(nbarValues.*StationaryDist.pdf))))/sum(sum(sum(sum(StationaryDist.pdf))));


%%
TFP_pdf=shiftdim(ValuesOnGrid(12,:,:,:,:),1);
ShareOfTFP(1)=nansum(nansum(TFP_pdf(logical((nbarValues>0).*(nbarValues<5))).*(StationaryDist.pdf(logical((nbarValues>0).*(nbarValues<5))))))/sum(sum(sum(sum(StationaryDist.pdf(logical((nbarValues>0).*(nbarValues<5)))))));
ShareOfTFP(2)=nansum(nansum(TFP_pdf(logical((nbarValues>=5).*(nbarValues<50))).*(StationaryDist.pdf(logical((nbarValues>=5).*(nbarValues<50))))))/sum(sum(sum(sum(StationaryDist.pdf(logical((nbarValues>=5).*(nbarValues<50)))))));
ShareOfTFP(3)=nansum(nansum(TFP_pdf(logical(nbarValues>=50)).*(StationaryDist.pdf(logical(nbarValues>=50)))))/sum(sum(sum(sum(StationaryDist.pdf(logical(nbarValues>=50))))));
ShareOfTFP(4)=nansum(nansum(nansum(nansum(TFP_pdf(logical(nbarValues>=0)).*StationaryDist.pdf(nbarValues>=0)))));


TFP_ear = sum(sum(sum(TFP_pdf(:,:,2,:).*(StationaryDist.pdf(:,:,2,:)))))/sum(sum(sum(StationaryDist.pdf(:,:,2,:))));

TFP_nonear = sum(sum(sum(TFP_pdf(:,:,1,:).*(StationaryDist.pdf(:,:,1,:)))))/sum(sum(sum(StationaryDist.pdf(:,:,1,:))));

%%
SUBStationaryDist.pdf=squeeze(StationaryDist.pdf(:,:,2,:));
nbarValues=squeeze(nbarValues(:,:,2,:));

SUBShareOfEstablishments(1)=100*(sum(sum(sum(SUBStationaryDist.pdf(logical((nbarValues>0).*(nbarValues<5))))))/sum(sum(sum(SUBStationaryDist.pdf))));
SUBShareOfEstablishments(2)=100*(sum(sum(sum(SUBStationaryDist.pdf(logical((nbarValues>=5).*(nbarValues<50))))))/sum(sum(sum(SUBStationaryDist.pdf))));
SUBShareOfEstablishments(3)=100*(sum(sum(sum(SUBStationaryDist.pdf(logical(nbarValues>=50)))))/sum(sum(sum(SUBStationaryDist.pdf))));
SUBShareOfEstablishments(4)=100*(sum(sum(sum(SUBStationaryDist.pdf))/sum(sum(sum(SUBStationaryDist.pdf)))));

SUBOutput_pdf=shiftdim(ProbDensityFns(2,:,:,2,:),1);
SUBShareOfOutput(1)=100*(sum(sum(sum(SUBOutput_pdf(logical((nbarValues>0).*(nbarValues<5))))))/sum(sum(sum(SUBOutput_pdf))));
SUBShareOfOutput(2)=100*(sum(sum(sum(SUBOutput_pdf(logical((nbarValues>=5).*(nbarValues<50))))))/sum(sum(sum(SUBOutput_pdf))));
SUBShareOfOutput(3)=100*(sum(sum(sum(SUBOutput_pdf(logical(nbarValues>=50))))))/sum(sum(sum(SUBOutput_pdf)));
SUBShareOfOutput(4)=100*(sum(sum(sum(SUBOutput_pdf)))/sum(sum(sum(SUBOutput_pdf))));

SUBLabour_pdf=shiftdim(ProbDensityFns(3,:,:,2,:),1);
SUBShareOfLabour(1)=100*(sum(sum(sum(SUBLabour_pdf(logical((nbarValues>0).*(nbarValues<5))))))/sum(sum(sum(SUBLabour_pdf))));
SUBShareOfLabour(2)=100*(sum(sum(sum(SUBLabour_pdf(logical((nbarValues>=5).*(nbarValues<50))))))/sum(sum(sum(SUBLabour_pdf))));
SUBShareOfLabour(3)=100*(sum(sum(sum(SUBLabour_pdf(logical(nbarValues>=50)))))/sum(sum(sum(SUBLabour_pdf))));
SUBShareOfLabour(4)=100*(sum(sum(sum(SUBLabour_pdf)))/sum(sum(sum(SUBLabour_pdf))));


SUBCapital_pdf=shiftdim(ProbDensityFns(1,:,:,2,:),1);
SUBShareOfCapital(1)=100*(sum(sum(sum(SUBCapital_pdf(logical((nbarValues>0).*(nbarValues<5))))))/sum(sum(sum(SUBCapital_pdf))));
SUBShareOfCapital(2)=100*(sum(sum(sum(SUBCapital_pdf(logical((nbarValues>=5).*(nbarValues<50))))))/sum(sum(sum(SUBCapital_pdf))));
SUBShareOfCapital(3)=100*(sum(sum(sum(SUBCapital_pdf(logical(nbarValues>=50)))))/sum(sum(sum(SUBCapital_pdf))));


%% Plots
plot_earmarkedmodel
