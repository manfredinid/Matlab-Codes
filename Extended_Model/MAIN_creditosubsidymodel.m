%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% Model presented in he TBrazilian Slump and the Government-driven 
% Credit Expansion (2020)

%% Initial setups
clear all;
close all;
Parallel=1; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
tic;

%% Model A
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=0.02; % Interest rate on earmarked credit
Params.g_ear=0.2; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.

creditsubsidymodel;

teste1 =(sum(sum(StationaryDist.pdf,3)));

% Agggregate Values
A.Output.Y=AggVars(2);
A.Output.N=AggVars(3);
A.Output.K=AggVars(1);
A.Output.KdivY=A.Output.K/A.Output.Y;
A.Output.TFP=(A.Output.Y/((A.Output.K^ Params.alpha)*(A.Output.N^ Params.gamma)));

% Agggregate Values without Subsidy
A.TAX.Output.Y=AggVars(8);
A.TAX.Output.N=AggVars(9);
A.TAX.Output.K=AggVars(7);
A.TAX.Output.KdivY=A.TAX.Output.K/A.TAX.Output.Y;
A.TAX.Output.TFP=(A.TAX.Output.Y/((A.TAX.Output.K^ Params.alpha)*(A.TAX.Output.N^ Params.gamma)));

% Agggregate Values with Subsidy
A.SUB.Output.Y=AggVars(5);
A.SUB.Output.N=AggVars(6);
A.SUB.Output.K=AggVars(4);
A.SUB.Output.KdivY=A.SUB.Output.K/A.SUB.Output.Y;
A.SUB.Output.TFP=(A.SUB.Output.Y/((A.SUB.Output.K^ Params.alpha)*(A.SUB.Output.N^ Params.gamma)));

A.ShareOfEstablishments=ShareOfEstablishments;
A.ShareOfOutput=ShareOfOutput;
A.ShareOfLabour=ShareOfLabour;
A.ShareOfCapital=ShareOfCapital;
A.AverageEmployment=AverageEmployment;


A.Percentage_tax=Percentage_tax;

%% Model B 
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=0.02; % Interest rate on earmarked credit
Params.g_ear=0.4; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.

creditsubsidymodel;

teste2 =(sum(sum(StationaryDist.pdf,3)));

% Agggregate Values
B.Output.Y=AggVars(2);
B.Output.N=AggVars(3);
B.Output.K=AggVars(1);
B.Output.KdivY=B.Output.K/B.Output.Y;
B.Output.TFP=(B.Output.Y/((B.Output.K^ Params.alpha)*(B.Output.N^ Params.gamma)));

% Agggregate Values without Subsidy
B.TAX.Output.Y=AggVars(8);
B.TAX.Output.N=AggVars(9);
B.TAX.Output.K=AggVars(7);
B.TAX.Output.KdivY=B.TAX.Output.K/B.TAX.Output.Y;
B.TAX.Output.TFP=(B.TAX.Output.Y/((B.TAX.Output.K^ Params.alpha)*(B.TAX.Output.N^ Params.gamma)));

% Agggregate Values with Subsidy
B.SUB.Output.Y=AggVars(5);
B.SUB.Output.N=AggVars(6);
B.SUB.Output.K=AggVars(4);
B.SUB.Output.KdivY=B.SUB.Output.K/B.SUB.Output.Y;
B.SUB.Output.TFP=(B.SUB.Output.Y/((B.SUB.Output.K^ Params.alpha)*(B.SUB.Output.N^ Params.gamma)));

B.ShareOfEstablishments=ShareOfEstablishments;
B.ShareOfOutput=ShareOfOutput;
B.ShareOfLabour=ShareOfLabour;
B.ShareOfCapital=ShareOfCapital;
B.AverageEmployment=AverageEmployment;


B.Percentage_tax=Percentage_tax;



%%

fprintf('         Market Rate     Subsidized Rate\n')
fprintf('Model A %9.2f  %12.2f  \n', [A.Percentage_tax])
fprintf('Model B %9.2f  %12.2f  \n', [B.Percentage_tax])


fprintf('                Model A      Model B\n');
fprintf('Total Output  %8.2f  %8.2f \n', [A.Output.Y B.Output.Y])
fprintf('Labor         %8.2f  %8.2f \n', [A.Output.N B.Output.N])
fprintf('Capital       %8.2f  %8.2f \n', [A.Output.K B.Output.K])
fprintf('TFP           %8.2f  %8.2f \n',[ A.Output.TFP B.Output.TFP])
fprintf('Percentage of firms with\n')



%%

fprintf('Model A  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                               <5     5 to 49     >=50    total\n');
fprintf('Share of establishments  %8.2f  %8.2f  %8.2f  %8.2f  \n', A.ShareOfEstablishments);
fprintf('Share of output          %8.2f  %8.2f  %8.2f  %8.2f\n', A.ShareOfOutput);
fprintf('Share of labor          %8.2f  %8.2f  %8.2f  %8.2f\n', A.ShareOfLabour);
fprintf('Share of capital         %8.2f  %8.2f  %8.2f  %8.2f\n', A.ShareOfCapital);
fprintf('Share of employment      %8.2f  %8.2f  %8.2f  %8.2f\n', A.AverageEmployment);


fprintf('Model B  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                               <5     5 to 49     >=50    total\n');
fprintf('Share of establishments  %8.2f  %8.2f  %8.2f  %8.2f  \n', B.ShareOfEstablishments);
fprintf('Share of output          %8.2f  %8.2f  %8.2f  %8.2f\n', B.ShareOfOutput);
fprintf('Share of labor          %8.2f  %8.2f  %8.2f  %8.2f\n', B.ShareOfLabour);
fprintf('Share of capital         %8.2f  %8.2f  %8.2f  %8.2f\n', B.ShareOfCapital);
fprintf('Share of employment      %8.2f  %8.2f  %8.2f  %8.2f\n', B.AverageEmployment);
%%

fprintf('Model A  \n');
fprintf('                     Total  with r_ear   with r_market\n');
fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[A.Output.TFP A.SUB.Output.TFP A.TAX.Output.TFP]);
fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [A.Output.Y A.SUB.Output.Y A.TAX.Output.Y]);
fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [A.Output.N A.SUB.Output.N A.TAX.Output.N]);
fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [A.Output.K A.SUB.Output.K A.TAX.Output.K]);



fprintf('Model B  \n');
fprintf('                     Total  with r_ear   with r_market\n');
fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[B.Output.TFP B.SUB.Output.TFP B.TAX.Output.TFP]);
fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [B.Output.Y B.SUB.Output.Y B.TAX.Output.Y]);
fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [B.Output.N B.SUB.Output.N B.TAX.Output.N]);
fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [B.Output.K B.SUB.Output.K B.TAX.Output.K]);
%%
figure;
plot(s_grid,teste1,'r' )
hold on;
plot(s_grid,teste2,'b' )

toc;