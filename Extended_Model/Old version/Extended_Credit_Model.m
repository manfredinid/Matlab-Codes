% Draft for the Credit Imbalance model

clear all;
close all;
Parallel=0; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU


rng('default') % For reproducibility
tic;
%% Some Toolkit options 

vfoptions.parallel=Parallel;
simoptions.parallel=Parallel;
heteroagentoptions.verbose=1;
simoptions.agententryandexit=1;
simoptions.endogenousexit=0;

%% Parameters

% Preferences
Params.beta=0.96; % Discount rate

% Production fn
Params.alpha=0.3;  % Capital share
Params.gamma=0.5; % alpha + gamma must be ~= 1
Params.delta=0.05; % Depreciation rate of physical capital
Params.cf=0; % Fixed cost of production

% Firm entry and exit
Params.Ne=1; % mass of new potential new entrant distribution.
Params.ce=1; % Fixed cost of entry (this is a normalization)
Params.lambda=0.1; % Probability of firm exit

% The actual 'distortionary policy rate'
Params.taurate=0; % This is the rate for the tax.
Params.subsidyrate=0; % This is the rate for the subsidy.
Params.gcost=0.5;

%% Exogenous processes; productivity and subsidies
n_s= 10; %Firm-specific Productivity level
n_sub = 20; %credit subsidy (must be an even number)


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
mu_a = 0.4;      % Mean (a).
sigma_a = 0.2;  % Standard deviation (a).
% Non-earmarked (b)
mu_b = 4;       % Mean (b).
sigma_b = 1;    % Standard deviation (b).
sz = [n_sub/states, 1];  % Size vector.
tau_grid = reshape([normrnd(mu_a, sigma_a, sz), normrnd(mu_b, sigma_b, sz)],[2*sz(1),1]);

% Transition matrix (considering that productivity and subsidy are independent) 
n_z=[n_s,length(tau_grid)];
z_grid=[s_grid; tau_grid];
% independent chains
pi_z=kron( pi_s,eye(prod(n_sub)))'; % transition matrix for the exogenous z variables
for ii = 1: length(pi_z)
A = round(sum(pi_z(:,ii)),5);
if A == 1
else
   error('transition matrix is wrong')
end
end
pi_z=pi_z';


%% Check endogenous, exogenous and decision variables

n_a=10;  % 401-- Hopenhayn & Rogerson (1993) use 250, this is intended as an approximation of continuous variable, so I use more.
%%%%% FIX FIX a_grid
%a_grid=[linspace(0,100,9),((logspace(0,pi,n_a-9-1)-1)/(pi-1))*(5000-9)+9,10^6]'; % One less point in standard grid
a_grid = sort(unifrnd(1,1000,n_a,1), 'descend');
n_d=0; % None.
d_grid=[];

disp('sizes')
disp('vector(s) of endogenous state variables')
disp(n_a)
disp('vector(s) of exogenous state variable')
disp(n_z)
disp('vector(s) of decision variabes')
disp(n_d)


%% Potential draws for the pair (s, psi)

sdist = 1:round(500/n_s,1):500;
smean = 1.5;
ssigma = 1;
cumsum_pistar_s = logncdf(sdist,smean,ssigma);
pistar_s=(cumsum_pistar_s-[0,cumsum_pistar_s(1:end-1)])';


taudist = 0:1/(n_sub-1):1;
cumsum_pistar_tau = betacdf(taudist,.5,.4);
pistar_tau =(cumsum_pistar_tau-[0,cumsum_pistar_tau(1:end-1)]);


%figure(1)
%subplot(1,2,1);
%plot(taudist,cumsum_pistar_tau)
%title('Potential draws for psi')
%subplot(1,2,2);
%plot(entrantsdist,cumsum_pistar_s)
%title('Potential draws for s')

%% Return Function
% The model extention would be an adjustment cost 

% Exit is exogenous - include as another 'DiscountFactorParamNames'
Params.oneminuslambda=1-Params.lambda; % This is now the conditional probability of survival.
% lambda is the average observed exit percentage between 2007--2017 
% (https://sidra.ibge.gov.br/Tabela/2718#resultado)
DiscountFactorParamNames={'beta','oneminuslambda'};

ReturnFn=@(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,taurate,subsidyrate, cf,gcost) RR2008p_ReturnFn(aprime_val, a_val,s_val, tau_val, p,r, alpha,gamma,taurate,subsidyrate, cf,gcost);
ReturnFnParamNames={'p','r','alpha','gamma','taurate','subsidyrate','cf','gcost'};

%% CHECK (erase)
% Incumbents exit in the beginning of the period
Params.i=1/Params.beta-1; % This is standard general eqm result in complete market models, comes from consumption euler eqn together with requirements of stationary eqm.
% The net return to capital in equilibrium will thus be
Params.r=Params.i+Params.delta; % That the gross return is just 1/beta-1 and equals i (that the gross return to capital equals the interest rate is a requirement of capital market clearance in model)
Params.p=1; 
Params.Ne=0.5;

vfoptions.parallel=Parallel;
if vfoptions.parallel==2
    V0=zeros([n_a,n_z,'gpuArray']);
else
    V0=zeros([n_a,n_z]);
end
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,d_grid,...
    a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, vfoptions);

figure;
surf(V(:,:,1))
%% Aspects of entry/exit
% Entry is endogenous and exit exogenous

% Both entry and exit matter for stationary distribution of agents
% Note: Because they are not a default part of agent simulation, you need 
% to pass the entry/exit aspects as part of simoptions.

EntryExitParamNames.DistOfNewAgents={'upsilon'};
% Probability of being in tau category

%%%%% FIX FIX teste
teste = repmat(1,[n_z n_a]);
%teste = eye(n_z(1), n_z(2),n_a);
Params.upsilon=pistar_s.*(pistar_tau).*teste;
Params.upsilon=Params.upsilon/sum(sum(sum(Params.upsilon)));

if (round(sum(sum(sum(Params.upsilon))),5) ~= 1)
    error('Upsilon is NOT a PMD.')
end

% Percentage of entering firms relative to existing agents
%Params.Ne=0.5;
EntryExitParamNames.MassOfNewAgents={'Ne'};

% Exogenous survival probability
EntryExitParamNames.CondlProbOfSurvival={'oneminuslambda'};

%% CHECK (erase)

% Check that everything is working so far by solving the simulation of agent distribution to get the stationary distribution.
simoptions.parallel=Parallel;
simoptions % Show which options are being set
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions,Params,EntryExitParamNames);

% Note: When using models, such as entry and exit, where the mass of agents is not equal to 1
% the toolkit will automatically keep track of distributions as StationaryDist.pdf and StationaryDist.mass

% Note: when tau=0, conditional on z there is a value of aprime that all
% firms will choose. This means that next period all firms will appear in
% just one of n_z values on a_grid (the pdf over on a_grid will appear as
% just n_z point masses; actually n_z+1 because of new entrants).
figure(2)
surf(StationaryDist.pdf(:,:,1))


%% Descriptions of SS values as functions
%Set initial values for prices
Params.p=1; 
Params.Ne=0.5;

GEPriceParamNames={'p','Ne'};

FnsToEvaluateParamNames(1).Names={'alpha','gamma'};
FnsToEvaluateFn_1 = @(aprime_val, s_val, tau_val, p,r, alpha,gamma,taurate)...
    (1-taurate)*p*s_val*(aprime_val^alpha)*(((((1-taurate)*s_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)))^gamma); % Total output
FnsToEvaluate={FnsToEvaluateFn_1};

%FnsToEvaluateParamNames(1).Names={};
%FnsToEvaluate={};

AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel, simoptions, EntryExitParamNames);

heteroagentoptions.specialgeneqmcondn={0,'entry'};

GeneralEqmEqn_1 = @(AggVars,p) 1/AggVars-p; 

GeneralEqmEqnParamNames(2).Names={'beta','ce'};
GeneralEqmEqn_Entry = @(EValueFn,p,beta,ce) beta*EValueFn-ce;


%% Equilibrium conditions

% 1 - Euler Equations
% 2 - Free Entry
% 3 - Market Clearing
% 4 - Stationary Distribution
% 5 - Optimal Production
% 6 - Entry/Exit Policies

% 1/6 Euler Equation
% Consumer's problem - complete markets solution
Params.i=1/Params.beta-1; % This is standard general eqm result in complete market models, comes from consumption euler eqn together with requirements of stationary eqm.
% The net return to capital in equilibrium will thus be
Params.r=Params.i+Params.delta; % That the gross return is just 1/beta-1 and equals i (that the gross return to capital equals the interest rate is a requirement of capital market clearance in model)


%2/6 Free entry and 3/5 Labor Market Clearing
GeneralEqmEqns={GeneralEqmEqn_1,GeneralEqmEqn_Entry};

%%
%Use the toolkit to find the equilibrium price index



if vfoptions.parallel==2
    V0=zeros([n_a,n_z],'gpuArray');
else
    V0=zeros([n_a,n_z]);
end
n_p=0;


disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,p_eqm_index, GeneralEqmCondition]=HeteroAgentStationaryEqm_Case1(V0, 0,...
    n_a, n_z, 0, pi_z, [], a_grid, z_grid, ReturnFn, FnsToEvaluate,...
    GeneralEqmEqns, Params, DiscountFactorParamNames, ReturnFnParamNames,...
    FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);

Params.p=p_eqm.p;

%%
%Now that we have the GE price, let's calculate the optimal deciosions

% 5/6 Optimal production 6/6 Entry/Exit Policies
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,[],a_grid,z_grid, pi_z, ReturnFn,...
    Params, DiscountFactorParamNames, ReturnFnParamNames, vfoptions);

% 4/6 Stationary Distribution 
StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z, simoptions, Params, EntryExitParamNames);

% Impose the labour market clearance, which involves calculating Ne. 
% Find mass of entry that clears the labor market.
%%%%%% CHANGE HR 1993
FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate) (p*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val)*z1_val)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluate={FnsToEvaluateFn_nbar};
AggValues=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);
InitialNe=Params.Ne;
Params.Ne=1/AggValues; % AggValues is presently equal to Nbar. This line is imposing/satisfying the labour market clearance condition.
StationaryDist.mass=StationaryDist.mass*(Params.Ne/InitialNe); % Take advantage of linearity of the stationary distribution in new entrants distribution.

%% Tables


FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate) (alpha/r)^((1-gamma)/(1-gamma-alpha)) *(gamma)^(gamma/(1-gamma-alpha)) *(z1_val*(1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*z2_val))^(1/(1-alpha-gamma));
FnsToEvaluateParamNames(2).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_nbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*p*z2_val)*z1_val)^(1/(1-alpha-gamma)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma)^((1-alpha)/(1-gamma-alpha)); % which evaluates to Nbar in the aggregate
FnsToEvaluateParamNames(3).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,mass, alpha,gamma,r,p,taurate,subsidyrate) ((1-((z2_val>=0)*taurate+(z2_val<0)*subsidyrate)*p*z2_val))^((alpha+gamma)/(1-gamma-alpha))*z1_val^(1/(1-gamma-alpha)) *(alpha/r)^(alpha/(1-gamma-alpha)) *(gamma)^(gamma/(1-gamma-alpha));
FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_nbar, FnsToEvaluateFn_output};

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1_Mass(StationaryDist.pdf,StationaryDist.mass, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames,EntryExitParamNames, n_d, n_a, n_z, [], a_grid, z_grid, Parallel,simoptions);

ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z, [], a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

% s_grid.^(1/(1-Params.gamma-Params.alpha))
nbarValues=shiftdim(ValuesOnGrid(2,:,:,:),1);
nbarValues=shiftdim(ValuesOnGrid(2,:,:,:),1);
normalize_employment=nbarValues(1,1,2); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=nbarValues./normalize_employment;

Partion1Indicator=logical(nbarValues<10);
Partion2Indicator=logical((nbarValues>=10).*(nbarValues<500));
Partion3Indicator=logical(nbarValues>=500);

% Check that the following is equal to prod(n_z), so 300
 if sum(sum(Partion1Indicator+Partion2Indicator+Partion3Indicator)) ~= prod(n_z)
     error('error')
 end
 
ShareOfEstablishments(1)=sum(sum(StationaryDist.pdf(Partion1Indicator)));
ShareOfEstablishments(2)=sum(sum(StationaryDist.pdf(Partion2Indicator)));
ShareOfEstablishments(3)=sum(sum(StationaryDist.pdf(Partion3Indicator)));

Output_pdf=shiftdim(ProbDensityFns(3,:,:,:),1);
ShareOfOutput(1)=sum(sum(sum(Output_pdf(Partion1Indicator))));
ShareOfOutput(2)=sum(sum(sum(Output_pdf(Partion2Indicator))));
ShareOfOutput(3)=sum(sum(sum(Output_pdf(Partion3Indicator))));

Labour_pdf=shiftdim(ProbDensityFns(2,:,:,:),1);
ShareOfLabour(1)=sum(sum(sum(Labour_pdf(Partion1Indicator))));
ShareOfLabour(2)=sum(sum(sum(Labour_pdf(Partion2Indicator))));
ShareOfLabour(3)=sum(sum(sum(Labour_pdf(Partion3Indicator))));

Capital_pdf=shiftdim(ProbDensityFns(1,:,:,:),1);
ShareOfCapital(1)=sum(sum(sum(Capital_pdf(Partion1Indicator))));
ShareOfCapital(2)=sum(sum(sum(Capital_pdf(Partion2Indicator))));
ShareOfCapital(3)=sum(sum(sum(Capital_pdf(Partion3Indicator))));

AverageEmployment(1)=sum(sum(nbarValues(Partion1Indicator).*StationaryDist.pdf(Partion1Indicator)))/sum(sum(StationaryDist.pdf(Partion1Indicator)));
AverageEmployment(2)=sum(sum(nbarValues(Partion2Indicator).*StationaryDist.pdf(Partion2Indicator)))/sum(sum(StationaryDist.pdf(Partion2Indicator)));
AverageEmployment(3)=sum(sum(nbarValues(Partion3Indicator).*StationaryDist.pdf(Partion3Indicator)))/sum(sum(StationaryDist.pdf(Partion3Indicator)));


disp('Share of establishments');
disp('       <10       10 to 490    >=500');
disp(ShareOfEstablishments);
disp('Share of output');
disp('       <10       10 to 490    >=500');
disp( ShareOfOutput);
disp('Share of labour');
disp('       <10       10 to 490    >=500');
disp(ShareOfLabour);
disp('Share of capital ');
disp('       <10       10 to 490    >=500');
disp(ShareOfCapital);
disp('Average employment');
disp('       <10       10 to 490    >=500');
disp(AverageEmployment);
