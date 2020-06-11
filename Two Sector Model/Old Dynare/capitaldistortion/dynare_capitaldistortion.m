%%
%clear all;
close all;

betaa =0.9;  
varpi = 0.5;
al =  1;
ah =3;
sigmaa = 1;
deltaa = 0.025; 
alphaa = 0.45;

% periods in the plot
n=100;

save parameterfile  betaa varpi al ah sigmaa deltaa alphaa


psi_h0=0;
psi_l0=0;

psi_hF_obser=0.1;
psi_lF_obser=-0.1;

psi_hF_simul=0.2;
psi_lF_simul=-0.2;

%%
%psi_h0=1.027/1.008; %2011 2.7 p.m.
%psi_l0=1; % 2011 0.8 p.m.

%psi_hF_obser=1.0356/1.008; %2013 2.36 p.m. % 2017 3.56
%psi_lF_obser= 1; %2013 0.6 p.m.  % 2017 0.79
 
%psi_hF_simul=1;
%psi_lF_simul=1;

%%
%dynare model_observed
dynare capitalpsi_obser

load('capitalpsi_obser_results.mat')
y_obser = y(1:n);
i_obser = i(1:n);
k_obser = k(1:n);
c_obser = c(1:n);
TFP_obser = tfp(1:n);
%kl_obser = kl(1:n);
%kh_obser = kh(1:n);

%% 
%dynare model_simul
dynare capitalpsi_simul

load('capitalpsi_simul_results.mat')
y_simul = y(1:n);
i_simul = i(1:n);
k_simul = k(1:n);
c_simul = c(1:n);
TFP_simul = tfp(1:n);
%kl_simul = kl(1:n);
%kh_simul = kh(1:n);


%% Graphs

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',0.8);
subplot(2,2,1);
plot(c_obser./c_obser(1));
hold on
plot(c_simul./c_simul(1),'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,2);
plot(i_obser./i_obser(1));
hold on
plot(i_simul./i_simul(1),'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,3);
plot(TFP_obser./TFP_obser(1));
hold on
plot(TFP_simul./TFP_simul(1),'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['TFP'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,4);
plot(y_obser./y_obser(1));
hold on
plot(y_simul./y_simul(1),'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')

snapnow

% capital-output ratio
figure;
plot(k_obser./y_obser);
hold on
plot(k_simul./y_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (observed)'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


%subplot(2,1,1);
%plot(kl_obser);
%hold on
%plot(kl_simul,'r');
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%title(['low-tecnhology capital'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
%subplot(2,1,2);
%plot(kh_obser);
%hold on
%plot(kh_simul,'r');
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%title(['high-technology capital'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
%snapnow

%% Selected Graphs

% In construction
%saveas(gcf,'kl_kh.eps') % save plot - place after the figure code



