dynare model_observed
load('model_observed_results.mat')
y_obser = y;
i_obser = i;
k_obser = k;
c_obser = c;

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


subplot(2,2,2);
plot(i_obser,'r-','LineWidth',1);
hold on
plot(i_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(k_obser,'r-','LineWidth',1);
hold on
plot(k_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(y_obser,'r-','LineWidth',1);
hold on
plot(y_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
suptitle(['observed policy']);
snapnow

% capital-output ratio
figure;
plot(k_obser./y_obser,'r-','LineWidth',1);
hold on
plot(k_simul./y_simul,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (observed)'],'FontSize',8,'FontWeight','bold');
snapnow




