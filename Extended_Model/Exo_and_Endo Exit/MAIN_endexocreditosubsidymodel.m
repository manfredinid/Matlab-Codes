%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% Model presented in the Brazilian Slump and the 
% Government-driven Credit Expansion (2020)


% higher exit cost higher number of small firms
%% Initial setups
clear all; clear mex; clear functions;clear java;
close all;
clearvars -global


Parallel=2; % 1 for (parallel) CPUs, 2 for GPU, 0 for single CPU
%vfoptions.lowmemory=1;
tic;

%% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

% Model A
A.Params.r_ear=(1+0.12)^(1/4)-1; % Interest rate on earmarked credit
A.Params.g_ear=0.4336;
A.Params.r_international= (1+0.2142)^(1/4)-1;


% Model B
B.Params.r_ear= (1+0.11)^(1/4)-1; % Interest rate on earmarked credit
B.Params.g_ear=0.5031;
B.Params.r_international =(1+0.2216)^(1/4)-1 ;

% Model C
C.Params.r_ear=(1+0.15)^(1/4)-1; % Interest rate on earmarked credit
C.Params.g_ear=0;
C.Params.r_international = (1+0.15)^(1/4)-1;

% Model D
D.Params.r_ear=(1+0.12)^(1/4)-1; % Interest rate on earmarked credit
D.Params.g_ear=0.4336;
D.Params.r_international = Params.r_hh;

%% Model A
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=A.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=A.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = A.Params.r_international;

%% Initial Guesses
Params.p=0.2939; % output price
Params.Ne=0.0601; % total mass of new entrants

%%
fprintf(2,'\nModel A  \n');
exitcreditsubsidymodel;
saveas(gcf,'initial','epsc')
%%
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
A.lambda=lambda;
%A.probenter=probenter;
A.ProbnbarValues=ProbnbarValues;



%% Model B 
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=B.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=B.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = B.Params.r_international;

%% Initial Guesses
%Params.p = 0.527;
%Params.Ne=0.51; % total mass of new entrants
%%

fprintf(2,'\nModel B  \n');
exitcreditsubsidymodel;
saveas(gcf,'observed','epsc')

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
B.lambda=lambda;
%B.probenter=probenter;
B.ProbnbarValues=ProbnbarValues;

%% Model C
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=C.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=C.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = C.Params.r_international;

%% Initial Guesses
%Params.p = 0.54;
%Params.Ne=0.37; % total mass of new entrants
%%

fprintf(2,'\nModel C  \n');
exitcreditsubsidymodel;
saveas(gcf,'alternative','epsc')

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

C.notexit=lambda;
%C.probenter=probenter;
C.ProbnbarValues=ProbnbarValues;

%% Model D
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=D.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=D.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = D.Params.r_international;

%% Initial Guesses
%Params.p = 0.54;
%Params.Ne=0.37; % total mass of new entrants
%%

fprintf(2,'\nModel C  \n');
exitcreditsubsidymodel;
saveas(gcf,'alternative','epsc')

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

D.notexit=lambda;
%D.probenter=probenter;
D.ProbnbarValues=ProbnbarValues;
%%

fprintf(' n_a  %8.3f \n', n_a);

fprintf(' n_s  %8.3f \n', n_s);

fprintf('                      Model A    Model B   Model C\n');
fprintf('r earmarked         %8.3f  %8.3f %8.3f\n',[ A.Params.r_ear B.Params.r_ear C.Params.r_ear])
fprintf('r Market            %8.3f  %8.3f %8.3f\n',[ A.Params.r_international B.Params.r_international C.Params.r_international])
fprintf('Subsidized Firms      %8.3f  %8.3f %8.3f\n',[ A.Params.g_ear B.Params.g_ear C.Params.g_ear])

fprintf('                      Model A    Model B   Model C\n');
fprintf('Total Output        %8.3f  %8.3f  %8.3f\n', [A.Output.Y B.Output.Y C.Output.Y])
fprintf('Labor               %8.3f  %8.3f  %8.3f\n', [A.Output.N B.Output.N C.Output.N])
fprintf('Capital             %8.3f  %8.3f  %8.3f\n', [A.Output.K B.Output.K C.Output.K])
fprintf('TFP                 %8.3f  %8.3f  %8.3f\n',[ A.Output.TFP B.Output.TFP C.Output.TFP])
fprintf('Net Foreign Assets  %8.3f  %8.3f  %8.3f\n',[ A.K_nfa B.K_nfa C.K_nfa])
fprintf('Subsidy Cost        %8.3f  %8.3f  %8.3f\n',[ A.cost B.cost C.cost])


%%

fprintf(2,'\nModel A  \n'); 
fprintf('Distribution statistics of benchmallrk economy  \n');
fprintf('                              <5     5 to 49   >=50   total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfCapital);
fprintf('Average Employment       %8.3f  %8.3f  %8.3f   %8.3f\n', A.AverageEmployment);


%%

fprintf(2,'\nModel B  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <5      5 to 49   >=50    total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfCapital);


%%


fprintf(2,'\nModel C  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <5      5 to 49   >=50     total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfCapital);

%%

fprintf(2,'\nModel B Subsidized\n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <5      5 to 49   >=50    total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', B.SUBShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', B.SUBShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', B.SUBShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', B.SUBShareOfCapital);
%%

fprintf(2,'\nModel A  \n');
fprintf('                  with r_ear  with r_market  Total\n');
fprintf('TFP               %8.3f    %8.3f    %8.3f    \n',[A.SUB.Output.TFP A.TAX.Output.TFP A.Output.TFP ]);
fprintf('Aggregate output  %8.3f    %8.3f    %8.3f \n', [ A.SUB.Output.Y A.TAX.Output.Y A.Output.Y]);
fprintf('Aggregate labor   %8.3f    %8.3f    %8.3f \n', [ A.SUB.Output.N A.TAX.Output.N A.Output.N]);
fprintf('Aggregate capital %8.3f    %8.3f    %8.3f \n ', [ A.SUB.Output.K A.TAX.Output.K A.Output.K]);



fprintf(2,'\nModel B \n');
fprintf('                  with r_ear  with r_market  Total\n');
fprintf('TFP               %8.3f    %8.3f    %8.3f    \n',[ B.SUB.Output.TFP B.TAX.Output.TFP B.Output.TFP]);
fprintf('Aggregate output  %8.3f    %8.3f    %8.3f \n', [ B.SUB.Output.Y B.TAX.Output.Y B.Output.Y]);
fprintf('Aggregate labor   %8.3f    %8.3f    %8.3f \n', [ B.SUB.Output.N B.TAX.Output.N B.Output.N]);
fprintf('Aggregate capital %8.3f    %8.3f    %8.3f \n ', [ B.SUB.Output.K B.TAX.Output.K B.Output.K]);

fprintf(2,'\nModel C \n');
fprintf('                  with r_ear  with r_market  Total\n');
fprintf('TFP               %8.3f    %8.3f    %8.3f    \n',[ C.SUB.Output.TFP C.TAX.Output.TFP C.Output.TFP]);
fprintf('Aggregate output  %8.3f    %8.3f    %8.3f \n', [ C.SUB.Output.Y C.TAX.Output.Y C.Output.Y]);
fprintf('Aggregate labor   %8.3f    %8.3f    %8.3f \n', [C.SUB.Output.N C.TAX.Output.N  C.Output.N]);
fprintf('Aggregate capital %8.3f    %8.3f    %8.3f \n ', [C.SUB.Output.K C.TAX.Output.K  C.Output.K]);
%%

fprintf('Average TFP  \n');
fprintf('                  <5         5 to 49      >=50       total \n');
fprintf('Model A       %8.3f     %8.3f    %8.3f   %8.3f \n',[A.ShareOfTFP])
fprintf('Model B       %8.3f     %8.3f    %8.3f   %8.3f\n',[B.ShareOfTFP])
fprintf('Model C       %8.3f     %8.3f    %8.3f   %8.3f\n',[C.ShareOfTFP])

%%

fprintf('Min TFP  \n');
fprintf('                  <5         5 to 49      >=50 \n');
fprintf('Model A       %8.3f     %8.3f    %8.3f    \n',[A.MinOfTFP])
fprintf('Model B       %8.3f     %8.3f    %8.3f    \n',[B.MinOfTFP])
fprintf('Model C       %8.3f     %8.3f    %8.3f    \n',[C.MinOfTFP])
%%
fprintf('                    Model A    Model B   Model C\n');
fprintf('Prob Stay         %8.3f  %8.3f %8.3f\n',[ A.notexit B.notexit C.notexit])
fprintf('Prob Enter        %8.3f  %8.3f %8.3f\n',[ A.probenter B.probenter C.probenter])
%%
figure;
%subplot(1,2,1)
%plot(s_grid,A.ProbnbarValues);
%xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,B.ProbnbarValues, ':k');
hold on;
plot(s_grid,C.ProbnbarValues', '-k');
xlim([0.9 2.5])
%hold on;
%line([s_grid(sum(C.ProbnbarValues==0)) s_grid(sum(C.ProbnbarValues==0))], get(gca, 'ylim'), 'Color', 'red',...
%    'LineStyle', '-');
%hold on;
%line([s_grid(sum(B.ProbnbarValues==0)) s_grid(sum(B.ProbnbarValues==0))], get(gca, 'ylim'), 'Color', 'black',...
%    'LineStyle', '-');
%title('earmarked')
xlabel('productivity')
%ylabel('employees')
legend('Final','Counterfactual', 'Location', 'northwest')

saveas(gcf,'proddist','epsc')

%

toc;