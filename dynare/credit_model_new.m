dynare model_observed
%% Graphs

figure;
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(c,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(i,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(k,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(y,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
suptitle(['observed policy']);
snapnow

% capital-output ratio
figure;
plot(k./y,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (observed)'],'FontSize',8,'FontWeight','bold');
snapnow


dynare model_simul
%% Graphs

figure;
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(c,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(i,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(k,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(y,'r-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
suptitle(['alternative policy']);
snapnow

% capital-output ratio
figure;
plot(k./y,'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (alternative)'],'FontSize',8,'FontWeight','bold');
snapnow


%load -mat modelsimul
%load -mat modelobser

