%% simple model based on Hopenhayn (1992)

vfoptions.parallel=2
simoptions.parallel=2

%% Parameters

% Initial guesses
Params.p=1; % output price
Params.Ne=0.8; % total mass of new entrants

% demand
Params.Dbar=400;
Params.beta=0.8;

% Firm-level technology
Params.alpha=0.3;  % Capital share
Params.cf=0; % Fixed cost of production% Parameters Calibration

% Entry and Exit
Params.ce=40; % Fixed cost of entry 
Params.lambda=0.1; % Probability of firm exit
Params.oneminuslambda=1-Params.lambda; % Probability of survival

% Declare discount factors
DiscountFactorParamNames={'beta','oneminuslambda'};
% Declare percentage of entrants
EntryExitParamNames.MassOfNewAgents={'Ne'};
% Exogenous survival probability
EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};

%% Exogenous state variables
% Exogenous AR(1) process on (log) productivity
% logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)
tauchenoptions.parallel=vfoptions.parallel;
Params.q=4;
Params.rho=0.93; % Hopenhayn & Rogerson (1993)
Params.sigma_logz=sqrt(0.53); % Hopenhayn & Rogerson (1993)
Params.sigma_epsilon=sqrt((1-Params.rho)*((Params.sigma_logz)^2));
Params.a=0.078; % Hopenhayn & Rogerson (1993) do not report, but Martin Flodén figures out the following (pg 5): http://martinfloden.net/files/macrolab.pdf
n_z=200; % I here call z, what Hopenhayn & Rogerson (1993) call s. The choice of n_z=20 follows them.
Params.q=4; % Hopenhayn & Rogerson (1993) do not report (based on Table 4 is seems something around q=4 is used, otherwise don't get values of z anywhere near as high as 27.3. (HR1993 have typo and call the column 'log(s)' when it should be 's') 
[z_grid, pi_z]=TauchenMethod(Params.a,Params.sigma_epsilon^2,Params.rho,n_z,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
z_grid=exp(z_grid);



logn = lognrnd(1,0.5,1,n_z);
cumsum_pistar_s = cumsum(logn./sum(logn));
pistar_z=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)])';
dist=1;
while dist>10^(-9)
    pistar_z_old=pistar_z;
    pistar_z=(pi_z')*pistar_z;
    dist=max(abs(pistar_z-pistar_z_old));
end

pi_z = eye(n_z,n_z);

%% Decision variables
%There is no d variable
n_a=1;
a_grid=1;
d_grid=[]; 
n_d=0;

%% Solve the Value Function
ReturnFn=@(aprime_val, a_val, s_val, p, alpha, cf) simple_ReturnFn(aprime_val, a_val, s_val, p, alpha, cf);
ReturnFnParamNames={'p', 'alpha', 'cf'}; %It is important that these are in same order as they appear in 'Hopenhayn1992_ReturnFn'

%% Aspects of the Endogenous entry
% Exit is exogenous with probability lambda
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;

EntryExitParamNames.DistOfNewAgents={'pistar_z'};
Params.pistar_z=pistar_z;

%% Use the toolkit to find the equilibrium price index

GEPriceParamNames={'p', 'Ne'}; 

heteroagentoptions.specialgeneqmcondn={0,'entry'};

FnsToEvaluateParamNames(1).Names={'alpha','p'};
FnsToEvaluateFn_output =@(aprime_val,a_val,z_val,mass,alpha,p)...
z_val*(z_val*p*alpha)^(alpha/(1-alpha));

FnsToEvaluate={FnsToEvaluateFn_output};


%%


GeneralEqmEqnParamNames(1).Names={'Dbar'};
GeneralEqmEqn_GoodsMarket = @(AggVars,p,Dbar) AggVars-Dbar/p(1);

GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,p,beta,ce) beta*EValueFn-ce; % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.

GeneralEqmEqns={GeneralEqmEqn_GoodsMarket,GeneralEqmEqn_Entry};

%%

heteroagentoptions.verbose=1;
n_p=0;

if vfoptions.parallel==2
    V0=zeros([n_a,n_z],'gpuArray');
else
    V0=zeros([n_a,n_z]);
end

disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);


%% Value Function, Policy and Firm Distribution in GE

disp('Calculating various equilibrium objects')
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
[V,Policy]=ValueFnIter_Case1(n_d,n_a,n_z,[],a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions, Params, EntryExitParamNames);

%%
figure;
plot(z_grid, StationaryDist.pdf)
hold on
plot(z_grid, pistar_z, '-r')

FnsToEvaluateParamNames(1).Names={'alpha','p'};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z_val,mass,alpha,p)(z_val*p*alpha)^(1/(1-alpha));

FnsToEvaluateParamNames(2).Names={'alpha','p'};
FnsToEvaluateFn_s = @(aprime_val,a_val,z_val,mass,alpha,p)z_val;

FnsToEvaluateParamNames(3).Names={'alpha','p'};
FnsToEvaluateFn_output =@(aprime_val,a_val,z_val,mass,alpha,p)...
z_val*(z_val*p*alpha)^(alpha/(1-alpha));

FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_s,FnsToEvaluateFn_output};
    
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
%%
%ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1_mass(StationaryDist, Policy, FnsToEvaluate, Params,...
%    FnsToEvaluateParamNames,EntryExitParamNames, n_d, n_a, n_z,...
%    [], a_grid, z_grid, simoptions.parallel,simoptions);

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1_Mass(StationaryDist.mass,...
    Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames,...
    EntryExitParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions);


ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate,...
    Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,...
    simoptions.parallel,simoptions,EntryExitParamNames);


%%
figure;
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
subplot(2,2,1);
plot(shiftdim(ValuesOnGrid(1,1,:),1),shiftdim(ProbDensityFns(1,1,:),1),'-','LineWidth',1);
legend('Capital','Location', 'Best');
subplot(2,2,2);
plot(squeeze(ValuesOnGrid(2,1,:)),squeeze(ProbDensityFns(2,1,:)) ,'-.','LineWidth',1);
legend('Productivity','Location', 'Best');
subplot(2,2,3);
plot(squeeze(ValuesOnGrid(3,1,:)),squeeze(ProbDensityFns(3,1,:)) ,':','LineWidth',1);
set(gca,'Fontsize',8);
legend('Output','Location', 'Best');
subplot(2,2,4);
plot(z_grid,(ones(size(z_grid))*AggVars(1))/StationaryDist.mass)
hold on
plot(z_grid,(ones(size(z_grid))*AggVars(2))/StationaryDist.mass)
hold on
plot(z_grid,(ones(size(z_grid))*AggVars(3))/StationaryDist.mass)
legend('Capital','Productivity','Output','Location', 'Best');
%%
% AggVars(1) is the same as StationaryDist.pdf.*StationaryDist.mass*shiftdim(ValuesOnGrid(1,1,:),1)

z_grid = linspace(0,max(squeeze(ProbDensityFns(2,1,:))),6)';
z_grid = linspace(0,0.02,6)';
figure;
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
subplot(3,1,1);
plot(shiftdim(ValuesOnGrid(1,1,:),1),StationaryDist.pdf,'-','LineWidth',1);
hold on;
plot((ones(size(z_grid))*AggVars(1))/StationaryDist.mass,z_grid,'-r')
legend('Capital Distribution','Aggregate Average','Location', 'Best');
subplot(3,1,2);
plot(squeeze(ValuesOnGrid(2,1,:)),StationaryDist.pdf ,'-.','LineWidth',1);
hold on;
plot((ones(size(z_grid))*AggVars(2))/StationaryDist.mass,z_grid,'-r')
legend('Productivity Distribution','Aggregate Average','Location', 'Best');
subplot(3,1,3);
plot(squeeze(ValuesOnGrid(3,1,:)),StationaryDist.pdf ,':','LineWidth',1);
hold on;
plot((ones(size(z_grid))*AggVars(3))/StationaryDist.mass,z_grid,'-r')
set(gca,'Fontsize',8);
legend('Output Distribution','Aggregate Average','Location', 'Best');
