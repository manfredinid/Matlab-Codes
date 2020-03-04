%%
close all;

betaa =0.9;  
varpi = 0.5;
al =  1;
ah =1.5;
sigmaa = 1.5;
deltaa = 0.03; 
alphaa = 0.7;

save parameterfile  betaa varpi al ah sigmaa deltaa alphaa


%pi_p0=1.76;
%pi_g0=1.47;

%pi_pF_obser=1.4;
%pi_gF_obser=0.7;

%pi_pF_simul=0.8;
%pi_gF_simul=0.8;

%%
pi_p0=1.027/1.008; %2011 2.7 p.m.
pi_g0=1; % 2011 0.8 p.m.

pi_pF_obser=1.0356/1.008; %2013 2.36 p.m. % 2017 3.56
pi_gF_obser= 1; %2013 0.6 p.m.  % 2017 0.79
 
pi_pF_simul=1;
pi_gF_simul=1;

%%
dynare model_observed
load('model_observed_results.mat')
y_obser = y;
i_obser = i;
k_obser = k;
c_obser = c;

%% 
dynare model_simul
load('model_simul_results.mat')
y_simul = y;
i_simul = i;
k_simul = k;
c_simul = c;


%% Graphs

figure;
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(c_obser,'r-','LineWidth',1);
hold on
plot(c_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,2);
plot(i_obser,'r-','LineWidth',1);
hold on
plot(i_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,3);
plot(k_obser,'r-','LineWidth',1);
hold on
plot(k_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,4);
plot(y_obser,'r-','LineWidth',1);
hold on
plot(y_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')

snapnow

% capital-output ratio
figure;
plot(k_obser./y_obser,'r-','LineWidth',1);
hold on
plot(k_simul./y_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (observed)'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')

snapnow




