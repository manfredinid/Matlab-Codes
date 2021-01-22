%%  Brazilian Slump and the Government-driven Credit Expansion

%% Dictionary
%Params.r_ear: Interest rate on earmarked credit (quarterly) [0,1]
%Params.g_ear: Share of entrants who have access to earmarket credit [0,1] 
%Params.r_international: Interest rate on non-earmarked credit (quarterly) [0,1]
%Params.g_tau: subsidized credit and poor credit acess correlation [0,1]
% FIX FIX 0 is uncorrelated
%g_tau=1
%Percentage tax and good credit     0.00
%Percentage tax and bad credit      0.64
%Percentage sub and good credit     0.00
%Percentage sub and bad credit      0.36

%Params.ctau: Interest rate on credit constrained firms (quarterly) [0,1]
% to measure correlation between 0 and 1 the value is ctau*gtau
% Params.g_tau: correlation between poor credit access and subsidy
   % if 0 all non-subsidized firms are credit constrain (poor credit
   % access)

%% Initial setups
clear all; clear mex; clear functions;clear java;
close all;
clearvars -global


Parallel=2; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
tic;

%% Scenarios

% Benchmark Model (Model A)
A.Params.r_ear=(1+0.12)^(1/4)-1; 
A.Params.g_ear=0.4336;
A.Params.r_international= (1+0.2142)^(1/4)-1;
A.Params.ctau=(1+0.07)^(1/4)-1;
A.Params.g_tau=0.32/A.Params.g_ear; % era 0.32

% Uncorrelated Distortions (Model B)
B.Params.r_ear= (1+0.11)^(1/4)-1; 
B.Params.g_ear=0.5031;
B.Params.r_international =(1+0.2216)^(1/4)-1 ;
B.Params.ctau=0.0210;
B.Params.g_tau=0;

% Correlated Distortions (Model C)
C.Params.r_ear=(1+0.11)^(1/4)-1; 
C.Params.g_ear=0.5031;
C.Params.r_international = (1+0.2216)^(1/4)-1 ;
C.Params.ctau=0.0210;
C.Params.g_tau=1;

% International Capital Flows (Model D)
D.Params.r_ear=(1+0.12)^(1/4)-1; 
D.Params.g_ear=0.4336;
%D.Params.r_international =1/Params.beta-1;
D.Params.ctau=0.0210;
D.Params.g_tau=0.43;

%% Benchmark Model (Model A)
fprintf(2,'\nBenchmark Model (Model A) \n');

Params.r_ear=A.Params.r_ear; 
Params.g_ear=A.Params.g_ear; 
Params.r_international = A.Params.r_international;
Params.g_tau=A.Params.g_tau;
Params.ctau=A.Params.ctau;

% Initial Guesses
Params.p=0.9776;%0.4331; % output price
Params.Ne=0.0138;%.0199; % total mass of new entrants
%%
% Model A
earmarkedmodel;
saveas(gcf,'modelA','epsc')

% Saved results

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
%A.MinOfTFP=MinOfTFP;

A.ShareOfEstablishments=ShareOfEstablishments;
A.ShareOfOutput=ShareOfOutput;
A.ShareOfLabour=ShareOfLabour;
A.ShareOfCapital=ShareOfCapital;
A.AverageEmployment=AverageEmployment;
A.ShareOfTFP=ShareOfTFP;

A.Percentage_tax=Percentage_tax;
A.K_nfa=K_nfa;
A.TFP_ear =AggVars(13);
A.TFP_nonear =AggVars(14);
A.ebar=Params.ebar;
A.lambda=ExitRateOfFirms;
%A.probenter=probenter;
%A.ProbnbarValues=ProbnbarValues;



fprintf(2,'\nModel A  \n'); 
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                              <5     5 to 49   >=50   total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfCapital);
fprintf('Average Employment       %8.3f  %8.3f  %8.3f   %8.3f\n', A.AverageEmployment);
fprintf('Exit Rate                %8.3f\n', A.lambda);

fprintf(2,'\n2010  \n'); 
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                              <5     5 to 49   >=50   total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', [41.03 51.59 7.38 100]);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', [1.35 8.45 90.2]);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', [3.35 26.42 70.24 100]);
%fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfCapital);
%fprintf('Average Employment       %8.3f  %8.3f  %8.3f   %8.3f\n', A.AverageEmployment);
fprintf('Exit Rate                %8.3f\n', (0.1625)^1/4);


%% Uncorrelated Distortions (Model B)

fprintf(2,'\nUncorrelated Distortions (Model B)\n');


Params.r_ear=B.Params.r_ear; 
Params.g_ear=B.Params.g_ear; 
Params.r_international = B.Params.r_international;
Params.g_tau=B.Params.g_tau;
Params.ctau=B.Params.ctau;

% Initial Guesses
%Params.p = 0.527;
%Params.Ne=0.51; % total mass of new entrants

% Model B
earmarkedmodel;
saveas(gcf,'modelB','epsc')

% Saved results

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
B.ShareOfTFP=ShareOfTFP;
%B.MinOfTFP=MinOfTFP;

B.SUBShareOfEstablishments=SUBShareOfEstablishments;
B.SUBShareOfOutput=SUBShareOfOutput;
B.SUBShareOfLabour=SUBShareOfLabour;
B.SUBShareOfCapital=SUBShareOfCapital;


B.Percentage_tax=Percentage_tax;
B.K_nfa=K_nfa;
B.TFP_ear =AggVars(13);
B.TFP_nonear =AggVars(14);
B.ebar=Params.ebar;
B.ExitRateOfFirms=ExitRateOfFirms;
%B.probenter=probenter;
%B.ProbnbarValues=ProbnbarValues;

%% Correlated Distortions (Model C)

fprintf(2,'\nCorrelated Distortions (Model C)\n');

Params.r_ear=C.Params.r_ear; 
Params.g_ear=C.Params.g_ear; 
Params.r_international = C.Params.r_international;
Params.g_tau=C.Params.g_tau;
Params.ctau=C.Params.ctau;

% Initial Guesses
%Params.p = 0.54;
%Params.Ne=0.37; % total mass of new entrants

% Model C
earmarkedmodel;
saveas(gcf,'modelC','epsc')

% Saved results

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
C.ShareOfTFP=ShareOfTFP;
%C.MinOfTFP=MinOfTFP;

C.SUBShareOfEstablishments=SUBShareOfEstablishments;
C.SUBShareOfOutput=SUBShareOfOutput;
C.SUBShareOfLabour=SUBShareOfLabour;
C.SUBShareOfCapital=SUBShareOfCapital;

C.Percentage_tax=Percentage_tax;
C.K_nfa=K_nfa;
C.TFP_ear =AggVars(13);
C.TFP_nonear =AggVars(14);
C.ebar=Params.ebar;

C.ExitRateOfFirms=ExitRateOfFirms;
%C.probenter=probenter;
%C.ProbnbarValues=ProbnbarValues;

%% International Capital Flows (Model D)

fprintf(2,'\nInternational Capital Flows (Model D)  \n');
D.Params.r_international =1/Params.beta-1;
earmarkedmodel;

Params.r_ear=D.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=D.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
%Params.r_international = D.Params.r_international;
Params.g_tau=D.Params.g_tau;
Params.ctau=D.Params.ctau;

% Initial Guesses
%Params.p = 0.54;
%Params.Ne=0.37; % total mass of new entrants

% Model D
earmarkedmodel;
saveas(gcf,'modelD','epsc')

% Saved results

% Agggregate Values
D.Output.Y=AggVars(2);
D.Output.N=AggVars(3);
D.Output.K=AggVars(1);
D.Output.KdivY=D.Output.K/D.Output.Y;
D.Output.TFP=(D.Output.Y/((D.Output.K^Params.alpha)*(D.Output.N^Params.gamma)));
D.cost = AggVars(11);

% Agggregate Values without Subsidy
D.TAX.Output.Y=AggVars(8);
D.TAX.Output.N=AggVars(9);
D.TAX.Output.K=AggVars(7);
D.TAX.Output.KdivY=D.TAX.Output.K/D.TAX.Output.Y;
D.TAX.Output.TFP=(D.TAX.Output.Y/((D.TAX.Output.K^ Params.alpha)*(D.TAX.Output.N^ Params.gamma)));

% Agggregate Values with Subsidy
D.SUB.Output.Y=AggVars(5);
D.SUB.Output.N=AggVars(6);
D.SUB.Output.K=AggVars(4);
D.SUB.Output.KdivY=D.SUB.Output.K/D.SUB.Output.Y;
D.SUB.Output.TFP=(D.SUB.Output.Y/((D.SUB.Output.K^ Params.alpha)*(D.SUB.Output.N^ Params.gamma)));

D.ShareOfEstablishments=ShareOfEstablishments;
D.ShareOfOutput=ShareOfOutput;
D.ShareOfLabour=ShareOfLabour;
D.ShareOfCapital=ShareOfCapital;
D.AverageEmployment=AverageEmployment;
D.ShareOfTFP=ShareOfTFP;
%D.MinOfTFP=MinOfTFP;

D.SUBShareOfEstablishments=SUBShareOfEstablishments;
D.SUBShareOfOutput=SUBShareOfOutput;
D.SUBShareOfLabour=SUBShareOfLabour;
D.SUBShareOfCapital=SUBShareOfCapital;

D.Percentage_tax=Percentage_tax;
D.K_nfa=K_nfa;
D.TFP_ear =AggVars(13);
D.TFP_nonear =AggVars(14);
D.ebar=Params.ebar;

D.ExitRateOfFirms=ExitRateOfFirms;
%D.probenter=probenter;
%D.ProbnbarValues=ProbnbarValues;
