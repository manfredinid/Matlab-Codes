% NEW graphs
% This program graph all the results of the model
t=1995:2017;
T=length(yt);


% 2) Normalized data t = 100 [t=7 (2001)]
for i = norm

At=exp(at);
teste = at - at(1);
At=At / At(i) * 100;
 
THt=exp(tht);
THt= THt / THt(i) * 100;
 
TKt=exp(tkt);
TKt=TKt / TKt(i) * 100;
 
TBt=exp(tbt);
TBt=TBt / TBt(i) * 100;
 

ytg=yt / yt(i) * 100;
ctg=ct/ ct(i) * 100;
htg=ht/ ht(i) * 100;
itg=it/ it(i) * 100;

end;

% 3) Plot Data on GDP, Consumption, Labor and Investiment
figure(1)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-')

subplot(2,2,1);
plot(t, ytg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
hold on;
vline(2010, 'r-');
xlim([1995 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');

subplot(2,2,2);
plot(t, ctg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
hold on;
vline(2010, 'r-');
xlim([1995 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');

subplot(2,2,3);
plot(t, htg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
hold on;
vline(2010, 'r-');
xlim([1995 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, itg ,'-','LineWidth',1);
hold on;
vline(2010, 'r-');
set(gca,'Fontsize',8);
xlim([1995 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle('Observed Data')



% 4) Plot all wedges (together)
figure(2);
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
plot(t, At ,'-','LineWidth',1);
hold on;
plot(t, THt ,'--','LineWidth',1);
hold on;
plot(t, TKt ,'-.','LineWidth',1);
hold on;
plot(t, TBt ,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2017]);
hold off;
legend('Productivity Wedge','Labor Wedge','Capital Wedge','Bond Wedge', 'Location', 'Best');
suptitle('Wedges')

indicator = [zeros(1,19) ones(1,3) zeros(1,1)]';
% 5) Plot Data on GDP, Consumption, Labor and Investiment (individual)
figure(3)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-')
shadedTimeSeries(t, ytg, indicator);
%plot(t, ytg ,'-','LineWidth',1);
%xlim([1995 2017]);
%ylim([80 140]);
title(['Output'],'FontSize',8,'FontWeight','bold');

figure(4)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-')
shadedTimeSeries(t, ctg, indicator);
%plot(t, ctg ,'-','LineWidth',1);
xlim([1995 2017]);
%ylim([80 140]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');

figure(5)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-')
shadedTimeSeries(t, htg, indicator);
%plot(t, htg ,'-','LineWidth',1);
xlim([1995 2017]);
%ylim([80 140]);
title(['Labor'],'FontSize',8,'FontWeight','bold');


figure(6)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-')
shadedTimeSeries(t, itg, indicator);
%plot(t, itg ,'-','LineWidth',1);
xlim([1995 2017]);
%ylim([80 140]);
title(['Investment'],'FontSize',8,'FontWeight','bold');



% 4) Plot all wedges (together)
figure;
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
plot(t(:,9:T), At(9:T,:) ,'-','LineWidth',1);
hold on;
plot(t(:,9:T), THt(9:T,:) ,'--','LineWidth',1);
hold on;
plot(t(:,9:T), TKt(9:T,:) ,'-.','LineWidth',1);
hold on;
plot(t(:,9:T), TBt(9:T,:) ,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
hold on;
vline(2014, 'r-');
hold off;
legend('Productivity Wedge','Labor Wedge','Capital Wedge','Bond Wedge', 'Location', 'Best');
suptitle('Wedges')

