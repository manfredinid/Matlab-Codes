dynare simple_model_kl_2
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
%suptitle(['all wedges']);
snapnow

% The k/y its not compatible with the values from the inventory capital 
% measure -- too high
figure;
plot(exp(k)./exp(y),'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio'],'FontSize',8,'FontWeight','bold');
snapnow
% The results are far away from the observed values

% the i/y is ok
figure;
plot(exp(i)./exp(y),'b-','LineWidth',1);
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment Output Ratio'],'FontSize',8,'FontWeight','bold');
snapnow

[c(2)-9.63 y(2)-9.87 i(2)-8.28]

[c(end)-9.60  y(end)-9.77 i(end)-7.9]



