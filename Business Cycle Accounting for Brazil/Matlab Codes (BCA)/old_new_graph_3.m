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

figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);

% 11) Plot Data on GDP, Consumption, Labor and Investiment without capital
% wedge
figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except bond wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_nob,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
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


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);

% 11) Plot Data on GDP, Consumption, Labor and Investiment without labor
% wedge
figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except labor wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_noh,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
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


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except capital wedge']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_nok,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
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
it_noa = ctfp_ex(2,:)/tfp_ex(2,i)*100;
end


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ytg,'-',t, yt_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['GDP'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, ctg,'-',t, ct_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, itg,'-',t, it_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);


figure;
set(gcf,'Color',[1,1,1]);
plot(t, htg,'-',t, ht_noa,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
title(['Labor'],'FontSize',8,'FontWeight','bold');
suptitle(['except tfp']);