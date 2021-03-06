%% FIRM DYNAMICS MODEL WITH CREDIT MISALLOCATION
% Model presented in the Brazilian Slump and the 
% Government-driven Credit Expansion (2020)

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
C.Params.r_ear=0; % Interest rate on earmarked credit
C.Params.g_ear=0;
C.Params.r_international = (1+0.15)^(1/4)-1;

%% Model A
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states

Params.r_ear=A.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=A.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = A.Params.r_international;
%%
Params.p=0.2939; % output price
Params.Ne=0.0694; % total mass of new entrants
%%
fprintf(2,'\nModel A  \n');
creditsubsidymodel;
saveas(gcf,'initial','epsc')

%% Agggregate Values
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
A.ShareOfTFP=ShareOfTFP;

A.Percentage_tax=Percentage_tax;
A.K_nfa=K_nfa;
A.TFP_ear =AggVars(13);
A.TFP_nonear =AggVars(14);
A.ebar=Params.ebar;

A.ProbnbarValues=ProbnbarValues;

A.NON_nbar = nanmean(nanmean(NONnbarValues(:,:,:),3));
A.SUB_nbar = nanmean(nanmean(SUBnbarValues(:,:,:),3));
%% Model B 
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states
%Params.p=0.2939; % output price
%Params.Ne=0.0694; % total mass of new entrants

Params.r_ear=B.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=B.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = B.Params.r_international;


fprintf(2,'\nModel B  \n');
creditsubsidymodel;
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

B.Percentage_tax=Percentage_tax;
B.K_nfa=K_nfa;
B.TFP_ear =AggVars(13);
B.TFP_nonear =AggVars(14);
B.ebar=Params.ebar;

B.ProbnbarValues=ProbnbarValues;
B.NON_nbar = nanmean(nanmean(NONnbarValues(:,:,:),3));
B.SUB_nbar = nanmean(nanmean(SUBnbarValues(:,:,:),3));
%% Model C
% Earmarked credit with embebed subsidies (psi)
% Exgoenous states


Params.r_ear=C.Params.r_ear; % Interest rate on earmarked credit
Params.g_ear=C.Params.g_ear; % Share of (unconditional) potential entrants who have access to earmarket credit. Note that conditional on entry this will not be same.
Params.r_international = C.Params.r_international;


fprintf(2,'\nModel C  \n');
creditsubsidymodel;
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

C.Percentage_tax=Percentage_tax;
C.K_nfa=K_nfa;
C.TFP_ear =AggVars(13);
C.TFP_nonear =AggVars(14);
C.ebar=Params.ebar;

C.ProbnbarValues=ProbnbarValues;

C.NON_nbar = nanmean(nanmean(NONnbarValues(:,:,:),3));
C.SUB_nbar = nanmean(nanmean(SUBnbarValues(:,:,:),3));
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
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                              <5     5 to 49   >=50   total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', A.ShareOfCapital);

A.ShareOfEstablishments' - [41.03 51.59 7.38 100]'
A.ShareOfLabour' - [3.35 26.42 70.24 100]'
%%

fprintf(2,'\nModel B  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <5      5 to 49   >=50    total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', B.ShareOfCapital);



fprintf(2,'\nModel C  \n');
fprintf('Distribution statistics of benchmark economy  \n');
fprintf('                             <5      5 to 49   >=50     total\n');
fprintf('Share of establishments  %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfEstablishments);
fprintf('Share of output          %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfOutput);
fprintf('Share of labor           %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfLabour);
fprintf('Share of capital         %8.3f  %8.3f  %8.3f   %8.3f\n', C.ShareOfCapital);


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


fprintf('Average TFP  \n');
fprintf('                  non-earmaked         earmakerd \n');
fprintf('Model A       %8.3f     %8.3f  \n',[A.TFP_nonear A.TFP_ear])
fprintf('Model B       %8.3f     %8.3f  \n',[B.TFP_nonear B.TFP_ear])
fprintf('Model C       %8.3f     %8.3f  \n',[C.TFP_nonear C.TFP_ear])


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
%xlim([0.9 2.5])
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

figure;
%subplot(1,2,1)
plot(s_grid,B.SUB_nbar,'-k');
%xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,B.NON_nbar,':k');
hold on;
plot(s_grid,C.NON_nbar,'--k');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('earmarked (ECE)','non-earmarked (ECE)', 'Counterfactual','Location', 'northwest')
%%
B = B.ebar;
B = double(B);

%%
figure; 
h=surf(B(:,:,1))
set(h,'LineStyle','none');
zticks([0 1])
xticks(1:n_s-1:n_s)
xticklabels({'1.4','2.75'});
yticks(1:n_a-1:n_a)
yticklabels({'0','1900000'});
set(h,'LineStyle','none');
xlabel('productivity')
ylabel('capital')
zlabel('entry decision, entry=1')
c = jet(2);
%colormap('gray');
colormap(c);
saveas(gcf,'2entryear','epsc')
figure; 
h=surf(B(:,:,2))
zticks([0 1])
xticks(1:n_s-1:n_s)
xticklabels({'1.4','2.75'});
yticks(1:n_a-1:n_a)
yticklabels({'0','1900000'});
set(h,'LineStyle','none');
xlabel('productivity')
ylabel('capital')
zlabel('entry decision, entry=1')
%colormap('gray');
colormap(c);
%colorbar
saveas(gcf,'2entrynonear','epsc')

%%

%%
figure; 
bar(s_grid,B(1,:,1))
%set(h,'LineStyle','none');
yticks([0 1])
%xticks(1:n_s-1:n_s)
%xticklabels({'1.5','2.75'});
%yticks(1:n_a-1:n_a)
%yticklabels({'0','1900000'});
%set(h,'LineStyle','none');
xlabel('productivity')
%ylabel('capital')
ylabel('entry decision, entry=1')
saveas(gcf,'entrynonear','epsc')
%%
figure; 
bar(s_grid,B(1,:,2))
yticks([0 1])
xlabel('productivity')
ylabel('entry decision, entry=1')
saveas(gcf,'entryear','epsc')

%%
toc;