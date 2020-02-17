%%
close all;

beta =0.92;  
varpi = 1/3;
al =  1;
ah =1.1;
sigma = 2;
delta = 0.03; 

save parameterfile  beta varpi al ah sigma delta


pi_p0=1.76;
pi_g0=1.47;

pi_pF_obser=1.4;
pi_gF_obser=0.7;

pi_pF_simul=0.8;
pi_gF_simul=0.8;

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




