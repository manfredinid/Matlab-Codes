
%% Find equilibrium prices

heteroagentoptions.verbose=1;
n_p=0;
% uncomment after erase the 'to be erase' chunks
% initial value function
%if vfoptions.parallel==2
%    V0=zeros([n_a,n_z,'gpuArray']);
%else
%    V0=zeros([n_a,n_z]);
%end

disp('Calculating price vector corresponding to the stationary eqm')
[p_eqm,p_eqm_index,GeneralEqmCondn]=HeteroAgentStationaryEqm_Case1(V0,...
    n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn,...
    FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames,...
    ReturnFnParamNames, FnsToEvaluateParamNames, GeneralEqmEqnParamNames,...
    GEPriceParamNames,heteroagentoptions, simoptions, vfoptions, EntryExitParamNames);
%% Value Function, Policy and Firm Distribution in GE

disp('Calculating various equilibrium objects')
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
[V,Policy]=ValueFnIter_Case1(V0, n_d,n_a,n_z,[],a_grid,z_grid, pi_z,...
    ReturnFn, Params, DiscountFactorParamNames, ReturnFnParamNames,vfoptions);

StationaryDist=StationaryDist_Case1(Policy,n_d,n_a,n_z,pi_z,...
    simoptions, Params, EntryExitParamNames);
%% Post GE values
FnsToEvaluateParamNames(1).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_kbar = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,...
    taurate,subsidyrate) aprime_val;
FnsToEvaluateParamNames(2).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_output = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,...
    r,p,taurate,subsidyrate) p*(1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*(aprime_val^alpha)*(...
    ((((1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)))^gamma);  

FnsToEvaluateParamNames(3).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_nbar =@(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate)...
 (((1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
 )*z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma));


FnsToEvaluateParamNames(4).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_subsidy =  @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,...
    r,p,taurate,subsidyrate) (z2_val<0)*p*(1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*(aprime_val^alpha)*(...
    ((((1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)))^gamma); 

FnsToEvaluateParamNames(5).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_tax =  @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,...
    r,p,taurate,subsidyrate) (z2_val>0)*p*(1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*(aprime_val^alpha)*(...
    ((((1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)))^gamma); 

FnsToEvaluateParamNames(6).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_none=  @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,...
    r,p,taurate,subsidyrate) (z2_val==0)*p*(1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*(aprime_val^alpha)*(...
    ((((1-(taurate*(z2_val>=0)-subsidyrate*(z2_val<0))...
    )*z1_val*p*gamma))^(1/(1-gamma)) *aprime_val^(alpha/(1-gamma)))^gamma); 

FnsToEvaluateParamNames(7).Names={'alpha','gamma','r','p','taurate','subsidyrate'};
FnsToEvaluateFn_taxsubnone=  @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,...
    r,p,taurate,subsidyrate) (z2_val<0)*2 + (z2_val==0)*3+(z2_val>0)*4; 

FnsToEvaluate={FnsToEvaluateFn_kbar, FnsToEvaluateFn_output,FnsToEvaluateFn_nbar,...
    FnsToEvaluateFn_subsidy,FnsToEvaluateFn_tax,FnsToEvaluateFn_none,FnsToEvaluateFn_taxsubnone};





%%
%FnsToEvaluateParamNames(1).Names={'alpha','gamma', 'delta','r','p','taurate','subsidyrate'};

%FnsToEvaluateParamNames(1).Names={};
% Capital
%FnsToEvaluateFn_capital = @(aprime_val,a_val,z1_val,z2_val,mass,alpha,gamma,r,p,taurate,subsidyrate) aprime_val; 
%FnsToEvaluate={FnsToEvaluateFn_capital};
  
AggVars=EvalFnOnAgentDist_AggVars_Case1(StationaryDist, Policy,...
    FnsToEvaluate, Params, FnsToEvaluateParamNames, n_d, n_a, n_z,...
    d_grid, a_grid, z_grid, simoptions.parallel,simoptions,EntryExitParamNames);

ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_Case1_Mass(StationaryDist.pdf,...
    StationaryDist.mass, Policy, FnsToEvaluate, Params,...
    FnsToEvaluateParamNames,EntryExitParamNames, n_d, n_a, n_z,...
    [], a_grid, z_grid, Parallel,simoptions);


ProbDensityFns=EvalFnOnAgentDist_pdf_Case1(StationaryDist, Policy, FnsToEvaluate,...
    Params, FnsToEvaluateParamNames, n_d, n_a, n_z, d_grid, a_grid, z_grid,...
    simoptions.parallel,simoptions,EntryExitParamNames);
%% Agggregate Values
Output.Y=AggVars(2);
Output.N=AggVars(3);
Output.K=AggVars(1);
Output.KdivY=Output.K/Output.Y;
%% Average values

Output.perN=AggVars(3)/StationaryDist.mass;
Output.perK=AggVars(1)/StationaryDist.mass;



%%

Establishments_sub=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(7,:,:,:),1)==4))));
Establishments100_tax=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(7,:,:,:),1)==2))));
Establishments_none=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(7,:,:,:),1)==3))));
Percentage_tax = [Establishments100_tax    Establishments_none    Establishments_sub] ;


Non_producing=100*sum(sum(sum(StationaryDist.pdf(shiftdim(ValuesOnGrid(3,:,:,:),1)<=0))));
Non_total=[Non_producing sum(Percentage_tax)];
%%
Output.TFP=(Output.Y/Output.N)./((Output.K/Output.N)^Params.alpha);



%%%%%%%%%%%%%%%
nbarValues=shiftdim(ValuesOnGrid(3,:,:,:),1);
normalize_employment=min(nonzeros(nbarValues)); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=nbarValues./normalize_employment;


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

fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                               <5     5 to 49     >=50    total\n');
fprintf('Share of establishments  %8.2f  %8.2f  %8.2f  %8.2f  \n', ShareOfEstablishments);
fprintf('Share of output          %8.2f  %8.2f  %8.2f  %8.2f\n', ShareOfOutput);
fprintf('Share of labor          %8.2f  %8.2f  %8.2f  %8.2f\n', ShareOfLabour);
fprintf('Share of capital         %8.2f  %8.2f  %8.2f  %8.2f\n', ShareOfCapital);
fprintf('Share of employment      %8.2f  %8.2f  %8.2f  %8.2f\n', AverageEmployment);

%% Display some output about the solution

fprintf('The equilibrium output price is p=%.4f \n', Params.p)
fprintf('The equilibrium value for the mass of entrants is Ne=%.4f \n', Params.Ne)

fprintf('Total Output is Y=%.4f \n', Output.Y)
fprintf('Labor is n=%.4f \n', Output.N)
fprintf('Capital is k=%.4f \n', Output.K)
fprintf('Total Factor Productivity is TFP=%.4f \n', Output.TFP)
fprintf('Total Subsided Output is Ysub=%.4f \n', AggVars(4))


fprintf('Taxed Firms      No tax or subsidy   Subsidized Firms\n')
fprintf('%9.2f  %12.2f  %19.2f  \n',Percentage_tax )

fprintf('Firms without production        Total(just for checking)  \n' )
fprintf('%9.2f  %25.2f  \n', Non_total)
