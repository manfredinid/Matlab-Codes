%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% Model presented in the Brazilian Slump and the 
% Government-driven Credit Expansion (2020)

%% Initial setups
clear all;
close all;
Parallel=1; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
tic;

%% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

% Model A
A.Params.r_ear=0.0517; % Interest rate on earmarked credit
A.Params.g_ear=0.429;
A.Params.r_international= 0.0287;


% Model B
B.Params.r_ear=0.0633; % Interest rate on earmarked credit
B.Params.g_ear=0.503;
B.Params.r_international = 0.0310;

% Model C
C.Params.r_ear=(1+0.13)^(1/4)-1; % Interest rate on earmarked credit
C.Params.g_ear=0;
C.Params.r_international = (1+0.13)^(1/4)-1;

%% Model A
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=A.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=A.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = A.Params.r_international;

creditsubsidymodel;

% Agggregate Values
A.Output.Y=AggVars(2);
A.Output.N=AggVars(3);
A.Output.K=AggVars(1);
A.Output.KdivY=A.Output.K/A.Output.Y;
A.Output.TFP=(A.Output.Y/((A.Output.K^Params.alpha)*(A.Output.N^Params.gamma)));
A.cost = AggVars(11);

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
A.K_nfa=K_nfa;

%% Model B 
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=B.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=B.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = B.Params.r_international;

creditsubsidymodel;



% Agggregate Values
B.Output.Y=AggVars(2);
B.Output.N=AggVars(3);
B.Output.K=AggVars(1);
B.Output.KdivY=B.Output.K/B.Output.Y;
B.Output.TFP=(B.Output.Y/((B.Output.K^ Params.alpha)*(B.Output.N^ Params.gamma)));
B.cost = AggVars(11);

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
B.SUB.Output.TFP=(B.SUB.Output.Y/((B.SUB.Output.K^Params.alpha)*(B.SUB.Output.N^Params.gamma)));


B.ShareOfEstablishments=ShareOfEstablishments;
B.ShareOfOutput=ShareOfOutput;
B.ShareOfLabour=ShareOfLabour;
B.ShareOfCapital=ShareOfCapital;
B.AverageEmployment=AverageEmployment;


B.Percentage_tax=Percentage_tax;
B.K_nfa=K_nfa;

%% Model A
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=C.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=C.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = C.Params.r_international;

creditsubsidymodel;


% Agggregate Values
C.Output.Y=AggVars(2);
C.Output.N=AggVars(3);
C.Output.K=AggVars(1);
C.Output.KdivY=C.Output.K/C.Output.Y;
C.Output.TFP=(C.Output.Y/((C.Output.K^Params.alpha)*(C.Output.N^Params.gamma)));
C.cost = AggVars(11);

% Agggregate Values without Subsidy
C.TAX.Output.Y=AggVars(8);
C.TAX.Output.N=AggVars(9);
C.TAX.Output.K=AggVars(7);
C.TAX.Output.KdivY=C.TAX.Output.K/C.TAX.Output.Y;
C.TAX.Output.TFP=(C.TAX.Output.Y/((C.TAX.Output.K^ Params.alpha)*(C.TAX.Output.N^ Params.gamma)));

% Agggregate Values with Subsidy
C.SUB.Output.Y=AggVars(5);
C.SUB.Output.N=AggVars(6);
C.SUB.Output.K=AggVars(4);
C.SUB.Output.KdivY=C.SUB.Output.K/C.SUB.Output.Y;
C.SUB.Output.TFP=(C.SUB.Output.Y/((C.SUB.Output.K^ Params.alpha)*(C.SUB.Output.N^ Params.gamma)));

C.ShareOfEstablishments=ShareOfEstablishments;
C.ShareOfOutput=ShareOfOutput;
C.ShareOfLabour=ShareOfLabour;
C.ShareOfCapital=ShareOfCapital;
C.AverageEmployment=AverageEmployment;


C.Percentage_tax=Percentage_tax;
C.K_nfa=K_nfa;


%%

fprintf('                      Model A      Model B\n');
fprintf('r earmarked         %8.2f  %8.2f \n',[ A.Params.r_ear B.Params.r_ear C.Params.r_ear])
fprintf('r Market            %8.2f  %8.2f \n',[ A.Params.r_international B.Params.r_international C.Params.r_international])
fprintf('Subzided Firms      %8.2f  %8.2f \n',[ A.Params.g_ear B.Params.g_ear C.Params.g_ear])

fprintf('                      Model A      Model B\n');
fprintf('Total Output        %8.2f  %8.2f \n', [A.Output.Y B.Output.Y C.Output.Y])
fprintf('Labor               %8.2f  %8.2f \n', [A.Output.N B.Output.N C.Output.N])
fprintf('Capital             %8.2f  %8.2f \n', [A.Output.K B.Output.K C.Output.K])
fprintf('TFP                 %8.2f  %8.2f \n',[ A.Output.TFP B.Output.TFP C.Output.TFP])
fprintf('Net Foreign Assets  %8.2f  %8.2f \n',[ A.K_nfa B.K_nfa C.K_nfa])
fprintf('Subsidy Cost        %8.2f  %8.2f \n',[ A.cost B.cost C.cost])


%%

fprintf(2,'\nModel A  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                              <10     10 to 49   >=50   total\n');
fprintf('Share of establishments  %8.2f  %8.2f  %8.2f   %8.2f\n', A.ShareOfEstablishments);
fprintf('Share of output          %8.2f  %8.2f  %8.2f   %8.2f\n', A.ShareOfOutput);
fprintf('Share of labor           %8.2f  %8.2f  %8.2f   %8.2f\n', A.ShareOfLabour);
fprintf('Share of capital         %8.2f  %8.2f  %8.2f   %8.2f\n', A.ShareOfCapital);
fprintf('Share of employment      %8.2f  %8.2f  %8.2f   %8.2f\n', A.AverageEmployment);


fprintf(2,'\nModel B  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <10      10 to 49   >=50    total\n');
fprintf('Share of establishments  %8.2f  %8.2f  %8.2f   %8.2f\n', B.ShareOfEstablishments);
fprintf('Share of output          %8.2f  %8.2f  %8.2f   %8.2f\n', B.ShareOfOutput);
fprintf('Share of labor           %8.2f  %8.2f  %8.2f   %8.2f\n', B.ShareOfLabour);
fprintf('Share of capital         %8.2f  %8.2f  %8.2f   %8.2f\n', B.ShareOfCapital);
fprintf('Share of employment      %8.2f  %8.2f  %8.2f   %8.2f\n', B.AverageEmployment);

fprintf(2,'\nModel C  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <10      10 to 49   >=50    total\n');
fprintf('Share of establishments  %8.2f  %8.2f  %8.2f   %8.2f\n', C.ShareOfEstablishments);
fprintf('Share of output          %8.2f  %8.2f  %8.2f   %8.2f\n', C.ShareOfOutput);
fprintf('Share of labor           %8.2f  %8.2f  %8.2f   %8.2f\n', C.ShareOfLabour);
fprintf('Share of capital         %8.2f  %8.2f  %8.2f   %8.2f\n', C.ShareOfCapital);
fprintf('Share of employment      %8.2f  %8.2f  %8.2f   %8.2f\n', C.AverageEmployment);
%%

fprintf(2,'\nModel A  \n');
fprintf('                     Total  with r_ear   with r_market\n');
fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[A.Output.TFP A.SUB.Output.TFP A.TAX.Output.TFP]);
fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [A.Output.Y A.SUB.Output.Y A.TAX.Output.Y]);
fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [A.Output.N A.SUB.Output.N A.TAX.Output.N]);
fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [A.Output.K A.SUB.Output.K A.TAX.Output.K]);



fprintf(2,'\nModel B \n');
fprintf('                     Total  with r_ear   with r_market\n');
fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[B.Output.TFP B.SUB.Output.TFP B.TAX.Output.TFP]);
fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [B.Output.Y B.SUB.Output.Y B.TAX.Output.Y]);
fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [B.Output.N B.SUB.Output.N B.TAX.Output.N]);
fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [B.Output.K B.SUB.Output.K B.TAX.Output.K]);

fprintf(2,'\nModel C \n');
fprintf('                     Total  with r_ear   with r_market\n');
fprintf('TFP %22.2f  %8.2f  %8.2f    \n',[C.Output.TFP C.SUB.Output.TFP C.TAX.Output.TFP]);
fprintf('Aggregate output  %8.2f  %8.2f  %8.2f \n', [C.Output.Y C.SUB.Output.Y C.TAX.Output.Y]);
fprintf('Aggregate labor   %8.2f  %8.2f  %8.2f \n', [C.Output.N C.SUB.Output.N C.TAX.Output.N]);
fprintf('Aggregate capital %8.2f  %8.2f  %8.2f \n ', [C.Output.K C.SUB.Output.K C.TAX.Output.K]);
%%
%figure;
%plot(s_grid,teste1,'r' )
%hold on;
%plot(s_grid,teste2,'b' )

% CAPITAL
%figure; plot((squeeze(sum(ProbDensityFns(1,:,:,:)))))

%OUTPUT
%figure; plot((squeeze(sum(ProbDensityFns(2,:,:,:)))))

%LABOR
%figure; plot((squeeze(sum(ProbDensityFns(3,:,:,:)))))



toc;