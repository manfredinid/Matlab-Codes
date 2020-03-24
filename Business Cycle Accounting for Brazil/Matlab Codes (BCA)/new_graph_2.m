% New graphs 2 -- only forescast

%xt = [zeros(3,T); at', tht'; tkt'; tbt'] all wedges
%xt = [zeros(3,T); at', zeros(3,T)] only TFP
%xt = [zeros(3,T); zeros(1,T); tht'; zeros(2,T)] only labor
%xt = [zeros(3,T); zeros(2,T); tkt'; zeros(1,T)] only capital
%xt = [zeros(3,T); zeros(3,T); tbt'] only bond
% Matrix xt is on the bca_wedge_alternative2 code (line 22)

% gx --- consumption, investiment, labor, gdp


% T is the yt length
T=length(yt);


% 1) TFP Wedge

xt = [zeros(3,T); at'; zeros(3,T)];
xt= [xt,zeros(7,1)];
TFP = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    TFP(:,j) = zt(:,j);
        
end

% 2) Labor Wedge

xt = [zeros(3,T); zeros(1,T); tht'; zeros(2,T)];
xt= [xt,zeros(7,1)];
labor = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    labor(:,j) = zt(:,j);
        
end

% 1) Capital Wedge

xt = [zeros(3,T); zeros(2,T); tkt'; zeros(1,T)];
xt= [xt,zeros(7,1)];
capital = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    capital(:,j) = zt(:,j);
        
end

% 4) Bond Wedge

xt = [zeros(3,T); zeros(3,T); tbt'];
xt= [xt,zeros(7,1)];
bond = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(:,j+1) = xtp(:);
       
    bond(:,j) = zt(:,j);
        
end

% 5) All Wedges

xt = [zeros(3,T); at'; tht'; tkt'; tbt'];
xt= [xt,zeros(7,1)];
all_wedges = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    all_wedges(:,j) = zt(:,j);
        
end

%6) Normalazed Wedges (forescast)

for i=norm;

TFP = exp(TFP);
TFP = TFP / TFP(i) * 100;

labor = exp(labor);
labor = labor / labor(i) * 100;

capital = exp(capital);
capital = capital / capital(i) * 100;

bond = exp(bond);
bond = bond / bond(i) * 100;

all_wedges = exp(all_wedges);
%all_wedges = (all_wedges / all_wedges(:,i)) * 100;
end

for i=norm
yt_all = all_wedges(4,:)/all_wedges(4,i)*100;
ht_all = all_wedges(3,:)/all_wedges(3,i)*100;
ct_all = all_wedges(1,:)/all_wedges(1,i)*100;
it_all = all_wedges(2,:)/all_wedges(2,i)*100;

yt_a = TFP(4,:)/TFP(4,i)*100;
ht_a = TFP(3,:)/TFP(3,i)*100;
ct_a = TFP(1,:)/TFP(1,i)*100;
it_a = TFP(2,:)/TFP(2,i)*100;

yt_h = labor(4,:)/labor(4,i)*100;
ht_h = labor(3,:)/labor(3,i)*100;
ct_h = labor(1,:)/labor(1,i)*100;
it_h = labor(2,:)/labor(2,i)*100;

yt_k = capital(4,:)/capital(4,i)*100;
ht_k = capital(3,:)/capital(3,i)*100;
ct_k = capital(1,:)/capital(1,i)*100;
it_k = capital(2,:)/capital(2,i)*100;

yt_b = bond(4,:)/bond(4,7)*100;
ht_b = bond(3,:)/bond(3,7)*100;
ct_b = bond(1,:)/bond(1,7)*100;
it_b = bond(2,:)/bond(2,7)*100;
end


t=t(:,9:T);
yt_all = yt_all(:,9:T);
ht_all = ht_all(:,9:T);
ct_all = ct_all(:,9:T);
it_all = it_all(:,9:T);

yt_a = yt_a(:,9:T);
ht_a = ht_a(:,9:T);
ct_a = ct_a(:,9:T);
it_a = it_a(:,9:T);

yt_h = yt_h(:,9:T);
ht_h = ht_h(:,9:T);
ct_h = ct_h(:,9:T);
it_h = it_h(:,9:T);

yt_k = yt_k(:,9:T);
ht_k = ht_k(:,9:T);
ct_k = ct_k(:,9:T);
it_k = it_k(:,9:T);

yt_b = yt_b(:,9:T);
ht_b = ht_b(:,9:T);
ct_b = ct_b(:,9:T);
it_b = it_b(:,9:T);

ytg = ytg(9:T,:);
htg = htg(9:T,:);
ctg = ctg(9:T,:);
itg = itg(9:T,:);

%% 7) Plot Data on GDP, Consumption, Labor and Investiment with all wedges
figure(8)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ctg,'-',t, ct_all,'r:','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, itg,'-',t, it_all,'r:','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, htg,'-',t, ht_all,'r:','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_all,'r:','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');

suptitle(['all wedges']);

%% 8) Plot Data on GDP, Consumption, Labor and Investiment with efficiency
% wedge
figure(9)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ctg,'-',t, ct_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, itg,'-',t, it_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, htg,'-',t, ht_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');

suptitle(['only efficiency wedge']);


%% 9) Plot Data on GDP, Consumption, Labor and Investiment with labor
% wedge
figure(10)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ctg,'-',t, ct_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, itg,'-',t, it_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, htg,'-',t, ht_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');

suptitle(['only labor wedge']);

%% 10) Plot Data on GDP, Consumption, Labor and Investiment with capital
% wedge
figure(11)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ctg,'-',t, ct_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, itg,'-',t, it_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, htg,'-',t, ht_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');

suptitle(['only capital wedge']);


%% 11) Plot Data on GDP, Consumption, Labor and Investiment with bond
% wedge
figure(12)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ctg,'-',t, ct_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
ylim([98 125])
title(['Consumption'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, itg,'-',t, it_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, htg,'-',t, ht_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');

suptitle(['only bond wedge']);


%% 12) Output and Wedges
%figure (13);
%%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(t, ytg,'-','LineWidth',1);
%hold on;
%plot(t, yt_a,'--','LineWidth',1);
%hold on;
%plot(t, yt_h,'-.','LineWidth',1);
%hold on;
%plot(t, yt_k,':','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold on;
%plot(t, yt_b,'-.o','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold off;
%legend('Output','Efficiency Wedge','Labor Wedge','Capital Wedge','Bond Wedge', 'Location', 'Best');
%suptitle('Otput an one Wedge')


%% 14) Consumption and Wedges
%figure (14);
%%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(t, ctg,'-','LineWidth',1);
%hold on;
%plot(t, ct_a,'--','LineWidth',1);
%hold on;
%plot(t, ct_h,'-.','LineWidth',1);
%hold on;
%plot(t, ct_k,':','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold on;
%plot(t, ct_b,'-.o','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold off;
%legend('Consumption','Efficiency Wedge','Labor Wedge','Capital Wedge','Bond Wedge', 'Location', 'Best');
%suptitle('Consumption and one Wedge')

%% 15) Labor and Wedges
%figure (15);
%%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(t, htg,'-','LineWidth',1);
%hold on;
%plot(t, ht_a,'--','LineWidth',1);
%hold on;
%plot(t, ht_h,'-.','LineWidth',1);
%hold on;
%plot(t, ht_k,':','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold on;
%plot(t, ht_b,'-.o','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold off;
%legend('Labor','Efficiency Wedge','Labor Wedge','Capital Wedge','Bond Wedge', 'Location', 'Best');
%suptitle('Labor and one Wedge')

%% 16) Labor and Wedges
%figure (16);
%%set(gcf,'Color',[1,1,1]);
%set(groot,'DefaultAxesColorOrder',[0 0 0],...
%      'DefaultAxesLineStyleOrder','-.|-|--|:')
%plot(t, itg,'-','LineWidth',1);
%hold on;
%plot(t, it_a,'--','LineWidth',1);
%hold on;
%plot(t, it_h,'-.','LineWidth',1);
%hold on;
%plot(t, it_k,':','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold on;
%plot(t, it_b,'-.o','LineWidth',1);
%set(gca,'Fontsize',8);
%xlim([2003 2017]);
%hold off;
%legend('Investment','Efficiency Wedge','Labor Wedge','Capital Wedge','Bond Wedge', 'Location', 'Best');
%suptitle('Investment and one Wedge')


%% 17) Plot Consumption
figure;
set(gcf,'Color',[1,1,1]);


subplot(2,2,1);
plot(t, ctg,'-',t, ct_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
ylim([98 128]);
legend('Consumption','Only efficiency wedge','Location', 'northwest');
%title(['Efficincy Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, ctg,'-',t, ct_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
ylim([98 128]);
legend('Consumption','Only capital wedge','Location', 'northwest')
%title(['Capital Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, ctg,'-',t, ct_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
ylim([98 128]);
legend('Consumption','Only labor wedge','Location', 'northwest')
%title(['Labor Wedge'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ctg,'-',t, ct_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
ylim([98 128]);
legend('Consumption','Only bond wedge','Location', 'northwest')
%title(['Bond Wedge'],'FontSize',8,'FontWeight','bold');

suptitle(['Consumption']);

%% 17) Plot Output
figure;
set(gcf,'Color',[1,1,1]);


subplot(2,2,1);
plot(t, ytg,'-',t, yt_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Output','Only efficiency wedge','Location', 'northwest');
%title(['Efficincy Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, ytg,'-',t, yt_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Output','Only capital wedge','Location', 'northwest')
%title(['Capital Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, ytg,'-',t, yt_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Output','Only labor wedge','Location', 'northwest')
%title(['Labor Wedge'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Output','Only bond wedge','Location', 'northwest')
%title(['Bond Wedge'],'FontSize',8,'FontWeight','bold');

suptitle(['Output']);

%% 17) Plot Labor
figure;
set(gcf,'Color',[1,1,1]);


subplot(2,2,1);
plot(t, htg,'-',t, ht_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Labor','Only efficiency wedge','Location', 'southwest');
%title(['Efficincy Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, htg,'-',t, ht_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Labor','Only capital wedge','Location', 'southwest')
%title(['Capital Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, htg,'-',t, ht_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Labor','Only labor wedge','Location', 'southwest')
%title(['Labor Wedge'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, htg,'-',t, ht_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Labor','Only bond wedge','Location', 'southwest')
%title(['Bond Wedge'],'FontSize',8,'FontWeight','bold');

suptitle(['Labor']);

%% 17) Plot Output
figure;
set(gcf,'Color',[1,1,1]);


subplot(2,2,1);
plot(t, itg,'-',t, it_a,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Investment','Only efficiency wedge','Location', 'northwest');
%title(['Efficincy Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,2);
plot(t, itg,'-',t, it_k,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Investment','Only capital wedge','Location', 'northwest')
%title(['Capital Wedge'],'FontSize',8,'FontWeight','bold');


subplot(2,2,3);
plot(t, itg,'-',t, it_h,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Investment','Only labor wedge','Location', 'northwest')
%title(['Labor Wedge'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, itg,'-',t, it_b,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%ylim([98 128]);
legend('Investment','Only bond wedge','Location', 'northwest')
%title(['Bond Wedge'],'FontSize',8,'FontWeight','bold');

suptitle(['Investment']);