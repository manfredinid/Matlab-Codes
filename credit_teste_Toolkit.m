% Draft for the Credit Imbalance model

clear all;
close all;
Parallel=0; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU


rng('default') % For reproducibility
tic;


%% Parameters

% Preferences
Params.beta=0.96; % Discount rate

% Production fn
Params.alpha=0.3;  % Capital share
Params.delta=0.05; % Depreciation rate of physical capital
Params.cf=0; % Fixed cost of production

% Firm entry and exit
Params.Ne=1; % mass of new potential new entrant distribution.
Params.ce=1; % Fixed cost of entry (this is a normalization)
Params.lambda=0.1; % Probability of firm exit


%% Exogenous processes; productivity and subsidies
n_s=10; %Firm-specific Productivity level
n_sub = 10; %credit subsidy (must be an even number)


% Exogenous AR(1) process on (log) productivity
% logz=a+rho*log(z)+epsilon, epsilon~N(0,sigma_epsilon^2)
Params.rho=0.93; 
Params.sigma_logz=sqrt(0.53); 
Params.sigma_epsilon=sqrt((1-Params.rho)*((Params.sigma_logz)^2));
Params.a=0.078; 

tauchenoptions.parallel=Parallel;
Params.q=4; % Hopenhayn & Rogerson (1993) do not report (based on Table 4 is seems something around q=4 is used, otherwise don't get values of z anywhere near as high as 27.3. (HR1993 have typo and call the column 'log(s)' when it should be 's') 
[s_grid, pi_s]=TauchenMethod(Params.a,Params.sigma_epsilon^2,Params.rho,n_s,Params.q,tauchenoptions); %[states, transmatrix]=TauchenMethod_Param(mew,sigmasq,rho,znum,q,Parallel,Verbose), transmatix is (z,zprime)
s_grid=exp(s_grid);

% Exogenous process on subsidy
% Random number from a bimodal distribution
states = 2;
% Earmarked (a)
mu_a = -1;      % Mean (a).
sigma_a = 0.5;  % Standard deviation (a).
% Non-earmarked (b)
mu_b = 2;       % Mean (b).
sigma_b = 1;    % Standard deviation (b).
sz = [n_sub/states, 1];  % Size vector.
x_2 = reshape([normrnd(mu_a, sigma_a, sz), normrnd(mu_b, sigma_b, sz)],[2*sz(1),1]);

figure(1)
subplot(1, 2, 1); hist(s_grid); title('Productivity');
subplot(1, 2, 2); histogram(x_2); title('Subsidy');

% Transition matrix (considering that productivity and subsidy are independent) 
n_z=[n_s,length(x_2)];
z_grid=[s_grid; x_2];
% independent chains
pi_z=kron( pi_s,eye(prod(n_sub)))'; % transition matrix for the exogenous z variables
for ii = 1: length(pi_z)
A = round(sum(pi_z(:,ii)),5);
if A == 1
else
   disp('transition matrix os wrong')
end
end
pi_z=pi_z';

%% Check endogenous, exogenous and decision variables

% Grids for a variables: there are none
n_a=1; % number of endogenous state variabels, if none na=1
a_grid=1;
n_d=0; 
d_grid=[];

disp('sizes')
disp('vector(s) of endogenous state variables')
disp(n_a)
disp('vector(s) of exogenous state variable')
disp(n_z)
disp('vector(s) of decision variabes')
disp(n_d)

%%%%%%%%% should I have this value for p or just the one for w 
Params.p=1; % "value of ce is chosen so that [Free entry condition] is satisfied with p=1, pg 930, Hopenhayn & Rogerson (1993)
Params.w=1; % Normalization
%% Distribution of potential entrants
%%%%%% CHANGE pistar_s is not sum to 1

entrantsdist = 1:round(500/n_sub,1):500;
entrantsmean = 0;
entrantssigma = 4;
cusum_pistar_z = logncdf(entrantsdist,entrantsmean,entrantssigma);
logncdf(10,entrantsmean,entrantssigma) %has to close to 0.78
cumsum_pistar_s = logncdf(entrantsdist,entrantsmean,entrantssigma);
pistar_s=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)])';
plot(entrantsdist,pistar_s)
%xlim([1 10])

%% Value Function Problem
% The model extention would be an adjustment cost 

% Exit is exogenous - include as another 'DiscountFactorParamNames'
Params.oneminuslambda=1-Params.lambda; % This is now the conditional probability of survival.
% lambda is the average observed exit percentage between 2007--2017 
% (https://sidra.ibge.gov.br/Tabela/2718#resultado)
DiscountFactorParamNames={'beta','oneminuslambda'};


% Incumbents exit in the beginning of the period
%%%%%%%% CHANGE has to adjust (mix of Hopenhayan and RR 2008 function)
%%%%% CHANGE ReturnFn AND ReturnFnParamNames
ReturnFn=@(n_val,aprime_val, a_val, s_val, p, alpha, cf) Hopenhayn1992_ReturnFn(n_val,aprime_val, a_val, s_val, p, alpha, cf);
ReturnFnParamNames={'p', 'alpha', 'cf'}; %It is important that these are in same order as they appear in 'Hopenhayn1992_ReturnFn'

% Check that everything is working so far by solving the value function
vfoptions.parallel=Parallel;
if vfoptions.parallel==2
    V0=zeros([n_a,n_z],'gpuArray');
else
    V0=zeros([n_a,n_z]);
end
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);


figure(3)
surf(shiftdim(V,1))
title('Value fn')

%% Stationary Distribution of Agents with entry and exit
% Entry is endogenous and exit exogenous
% Both entry and exit matter for stationary distribution of agents
% Note: Because they are not a default part of agent simulation, you need 
% to pass the entry/exit aspects as part of simoptions.

simoptions.agententryandexit=1;

% Aspects of entry/exit

% upsilon has to be a PMF
pistar_tau = unidpdf(1:n_sub,n_sub);
EntryExitParamNames.DistOfNewAgents={'upsilon'};
% Probability of being in tau category
Params.upsilon=pistar_s.*(pistar_tau);


% Percentage of entering firms relative to existing agents
%%%%%%%%% ???????????
% The initial guess for the mass of existing agents is always 1
Params.Ne=0.5;
EntryExitParamNames.MassOfNewAgents={'Ne'};

% Exogenous survival probability
EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};


simoptions.parallel=Parallel;
simoptions % Show which options are being set
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);


% If you wanted to look at the pdf:
figure;
 surf(shiftdim(StationaryDist.pdf,1))
%% General equilibrium conditions 

%%%%%%% CHECK CHECK
%Use the toolkit to find the equilibrium prices
GEPriceParamNames={'Ne'};

FnsToEvaluateParamNames(1).Names={};
% Note: With entry-exit the mass of the distribution of agents often
% matters. So it becomes an extra input arguement in all functions to be evaluated.
% FnsToEvaluateFn_1 = @(aprime_val,a_val,z_val,agentmass,alpha,p) p*z_val*(aprime_val^alpha); % Total output
FnsToEvaluate={};

% Just to test: (note, is same command as usual, just need to include the optional extra inputs 'simoptions' and 'EntryExitParamNames' which contains all the needed info about entry/exit)
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

% The general equilibrium condition is that the EV^e-ce=0.
% This does not fit standard format for general equilibrium conditions.
heteroagentoptions.specialgeneqmcondn={0,'entry'};
% Certain kinds of general equilibrium conditions that are non-standard can
% be used via heteroagentoptions.specialgeneqmcondn

GeneralEqmEqnParamNames(2).Names={'p'};
GeneralEqmEqn_Entry = @(EValueFn,GEprices,p) EValueFn-p*GEprices(1); % Free entry conditions (expected returns equal zero in eqm); note that the first 'General eqm price' is ce, the fixed-cost of entry.

% The entry condition looks slightly different to more standard @(EValueFn,p,params)
% This is because 'p' is the name of a parameter, and so have used 'GEprices'
% instead of my usual 'p' to refer to the general equilibrium prices (here 'ce' and 'Ne')
GeneralEqmEqns={GeneralEqmEqn_Entry};
% Note that GeneralEqmEqn_Entry needed to be pointed out as special because
% it depends on the distribution of entrants and not the distribution of
% existing agents (all standard general eqm conditions involve the later).

teste=0000

%%
heteroagentoptions.verbose=1;
n_p=0;
disp('Calculating price vector corresponding to the stationary eqm')
% tic;
% NOTE: EntryExitParamNames has to be passed as an additional input compared to the standard case.
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(V0, n_d, n_a, ...
    n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns,...
    Params, DiscountFactorParamNames, ReturnFnParamNames, FnsToEvaluateParamNames,...
    GeneralEqmEqnParamNames, GEPriceParamNames,heteroagentoptions, simoptions,...
    vfoptions, EntryExitParamNames);
% findeqmtime=toc

Params.ce=p_eqm.ce;
Params.Ne=p_eqm.Ne;
%%


teste =555555
toc;
% Economy constrain include net exports
% S+ NX = I
