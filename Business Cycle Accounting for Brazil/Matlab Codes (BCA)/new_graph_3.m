% New graphs 3 -- except forescast

%xt = [zeros(3,T); at'; tht'; tkt'; tbt'] all wedges

% gx --- consumption, investiment, labor, gdp

% T is the yt length
T=length(yt);


% 1) Except Bond Wedge

xt = [zeros(3,T); at'; tht'; tkt'; zeros(1,T)];
xt= [xt,zeros(7,1)];
bond_ex = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    bond_ex(:,j) = zt(:,j);
        
end

bond_ex = exp(bond_ex);
for i=norm
bond_ex = bond_ex / bond_ex(i) * 100;
yt_nob = bond_ex(4,:)/bond_ex(4,i)*100;
ht_nob = bond_ex(3,:)/bond_ex(3,i)*100;
ct_nob = bond_ex(1,:)/bond_ex(1,i)*100;
it_nob = bond_ex(2,:)/bond_ex(2,i)*100;
end

bond_ex = bond_ex(:,9:T);
yt_nob = yt_nob(:,9:T);
ht_nob = ht_nob(:,9:T);
ct_nob = ct_nob(:,9:T);
it_nob = it_nob(:,9:T);

figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);

% 11) Plot Data on GDP, Consumption, Labor and Investiment without capital
% wedge
figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);


% 2) Except Labor Wedge

xt = [zeros(3,T); at'; zeros(1,T); tkt'; tbt'];
xt= [xt,zeros(7,1)];
labor_ex = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    labor_ex(:,j) = zt(:,j);
        
end

labor_ex = exp(labor_ex);
for i=norm
labor_ex = labor_ex / labor_ex(i) * 100;
yt_noh = labor_ex(4,:)/labor_ex(4,i)*100;
ht_noh = labor_ex(3,:)/labor_ex(3,i)*100;
ct_noh = labor_ex(1,:)/labor_ex(1,i)*100;
it_noh = labor_ex(2,:)/labor_ex(2,i)*100;
end

labor_ex = labor_ex(:,9:T);
yt_noh = yt_noh(:,9:T);
ht_noh = ht_noh(:,9:T);
ct_noh = ct_noh(:,9:T);
it_noh = it_noh(:,9:T);

figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);

% 11) Plot Data on GDP, Consumption, Labor and Investiment without labor
% wedge
figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);


% 2) Except capital Wedge

xt = [zeros(3,T); at'; tht'; zeros(1,T); tbt'];
xt= [xt,zeros(7,1)];
capital_ex = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    capital_ex(:,j) = zt(:,j);
        
end

capital_ex = exp(capital_ex);
for i=norm
capital_ex = capital_ex / capital_ex(i) * 100;
yt_nok = capital_ex(4,:)/capital_ex(4,i)*100;
ht_nok = capital_ex(3,:)/capital_ex(3,i)*100;
ct_nok = capital_ex(1,:)/capital_ex(1,i)*100;
it_nok = capital_ex(2,:)/capital_ex(2,i)*100;
end

capital_ex = capital_ex(:,9:T);
yt_nok = yt_nok(:,9:T);
ht_nok = ht_nok(:,9:T);
ct_nok = ct_nok(:,9:T);
it_nok = it_nok(:,9:T);

figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);




% 2) Except efficiency Wedge

xt = [zeros(3,T); zeros(1,T); tht'; tkt'; tbt'];
xt= [xt,zeros(7,1)];
tfp_ex = zeros(4,T); 
zt = zeros(4,T);

for j=1:T
    
    zt(:,j) = gx * xt(:,j);
    
    xtp = hx * xt(:,j);
    
    xt(1:3,j+1) = xtp(1:3);
       
    tfp_ex(:,j) = zt(:,j);
        
end

tfp_ex = exp(tfp_ex);
for i=norm
tfp_ex = tfp_ex / tfp_ex(i) * 100;
yt_noa = tfp_ex(4,:)/tfp_ex(4,i)*100;
ht_noa = tfp_ex(3,:)/tfp_ex(3,i)*100;
ct_noa = tfp_ex(1,:)/tfp_ex(1,i)*100;
it_noa = tfp_ex(2,:)/tfp_ex(2,i)*100;
end


tfp_ex = tfp_ex(:,9:T);
yt_noa = yt_noa(:,9:T);
ht_noa = ht_noa(:,9:T);
ct_noa = ct_noa(:,9:T);
it_noa = it_noa(:,9:T);

figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
title(['Labor'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);









%% 12) Output and Wedges
figure (13);
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
plot(t, ytg,'-','LineWidth',1);
hold on;
plot(t, yt_noa,'--','LineWidth',1);
hold on;
plot(t, yt_noh,'-.','LineWidth',1);
hold on;
plot(t, yt_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%hold on;
%plot(t, yt_nob,'-.o','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
hold off;
legend('Output','no Efficiency Wedge','no Labor Wedge','no Capital Wedge', 'Location',  'northwest');
suptitle('Output an one Wedge')


%% 14) Consumption and Wedges
figure (14);
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
plot(t, ctg,'-','LineWidth',1);
hold on;
plot(t, ct_noa,'--','LineWidth',1);
hold on;
plot(t, ct_noh,'-.','LineWidth',1);
hold on;
plot(t, ct_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%hold on;
%plot(t, ct_nob,'-.o','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
hold off;
legend('Consumption','no Efficiency Wedge','no Labor Wedge','no Capital Wedge', 'Location', 'northwest');
suptitle('Consumption and one Wedge')

%% 15) Labor and Wedges
figure (15);
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
plot(t, htg,'-','LineWidth',1);
hold on;
plot(t, ht_noa,'--','LineWidth',1);
hold on;
plot(t, ht_noh,'-.','LineWidth',1);
hold on;
plot(t, ht_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%hold on;
%plot(t, ht_nob,'-.o','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
hold off;
legend('Labor','no Efficiency Wedge','no Labor Wedge','no Capital Wedge', 'Location', 'northwest');
suptitle('Labor and one Wedge')

%% 16) Labor and Wedges
figure;
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
plot(t, itg,'-','LineWidth',1);
hold on;
plot(t, it_noa,'--','LineWidth',1);
hold on;
plot(t, it_noh,'-.','LineWidth',1);
hold on;
plot(t, it_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
%hold on;
%plot(t, it_nob,'-.o','LineWidth',1);
set(gca,'Fontsize',8);
xlim([2003 2017]);
hold off;
legend('Investment','no Efficiency Wedge','no Labor Wedge','no Capital Wedge', 'Location', 'northwest');
suptitle('Investment and one Wedge')
