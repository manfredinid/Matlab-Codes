% Plot TVP - second 

clear all;
close all;
clc;
sensitivity = 0;

% Escolha dos anos da FIR
FIR1 = 1963;
FIR2 = 1979;
FIR3 = 2015;


% Hiperparameters k
% Set some hyperparameters here (see page Koop 831, end of section 4.1)
k_Q = 0.1;
k_S = 0.1;
k_W = 0.01;
%%
TVP_VAR_difusa;

if istore ==1
     
    qus = [.10, .5, .90];
    imp75XY=squeeze(quantile(imp75,qus));
    imp81XY=squeeze(quantile(imp81,qus));
    imp96XY=squeeze(quantile(imp96,qus));
end

save('Model1')

if sensitivity == 1

k_Q = 0.1;
k_S = 0.01;
k_W = 0.001;

if istore ==1
     
    qus = [.10, .5, .90];
    imp75XY=squeeze(quantile(imp75,qus));
    imp81XY=squeeze(quantile(imp81,qus));
    imp96XY=squeeze(quantile(imp96,qus));
end

    TVP_VAR_difusa;
    save('Model2')
end



load yearlab.dat
%% Figura das S�ries Originais
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
  
subplot(2,2,1)
plot(yearlab,Y(:,1),'-','LineWidth',0.7);
set(gca,'Fontsize',8);
title('GDP per capita (annual %)')
xlim([1954 2017])

subplot(2,2,2)
plot(yearlab,Y(:,2),'-','LineWidth',0.7);
set(gca,'Fontsize',8);
title('TFP')
xlim([1954 2017])

subplot(2,2,[3,4])
indicator = [ones(1,sum(yearlab<1964)) zeros(1,sum((yearlab>=1964).*(yearlab<1974)))...
    ones(1,sum((yearlab>=1974).*(yearlab<1980)))  zeros(1,sum((yearlab>=1980).*(yearlab<2009)))...
    ones(1,sum((yearlab>=2009).*(yearlab<2015)))  zeros(1,sum((yearlab>=2015).*(yearlab<=2017)))]';
shadedTimeSeries(yearlab, Y(:,3) ,[indicator], '', {''}, [0 0 0]+0.8, 07);
txt = {'\leftarrow Authoritarian national','developmentalism'};
text(1976, 15 ,txt,'FontSize',10)
txt = {'New','state-led development \rightarrow '};
text(2003, 5,txt,'FontSize',10,'HorizontalAlignment', 'center')
txt = {'\leftarrow Golden age of','import substitution'};
text(1960, 20 ,txt,'FontSize',10)
title('BNDES disbursements to GFCF')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 20]);
saveas(gcf,'allvars','epsc')

%% C�lculo das Fun��es Impulso Resposta
% Os anos das Fun��es est�o no script FIR

if istore == 1
    %%
    load('Model1')
    %%
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':|-|--|:')
 plot(1:nhor,squeeze(imp75XY(2,1,:)),'LineWidth',1.2)
    hold on
    plot(1:nhor,squeeze(imp81XY(2,1,:)),'LineWidth',1.2)
    hold on
    plot(1:nhor,squeeze(imp96XY(2,1,:)),'LineWidth',1.2)
    lgd=legend('1963','1979','2015');
    lgd.FontSize = 12;
set(gca,'Fontsize',8);
 xlim([1 nhor]);
%title(['Response Growth'],'FontSize',8,'FontWeight','bold');
saveas(gcf,'ResponseGrowth','epsc')


figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':|-|--|:')
 plot(1:nhor,squeeze(imp75XY(2,2,:)),'LineWidth',1.2)
    hold on
    plot(1:nhor,squeeze(imp81XY(2,2,:)),'LineWidth',1.2)
    hold on
    plot(1:nhor,squeeze(imp96XY(2,2,:)),'LineWidth',1.2)
 xlim([1 nhor]);
%title(['Response TFP'],'FontSize',8,'FontWeight','bold');
saveas(gcf,'ResponseTFP1','epsc')

figure;
  plot(1:nhor,squeeze(imp75XY(2,3,:)),'LineWidth',1.2)
    hold on
    plot(1:nhor,squeeze(imp81XY(2,3,:)),'LineWidth',1.2)
    hold on
    plot(1:nhor,squeeze(imp96XY(2,3,:)),'LineWidth',1.2)
 xlim([1 nhor]);
%title(['Response BNDES'],'FontSize',10,'FontWeight','bold');
saveas(gcf,'ResponseBNDES','epsc')
end
%%
load yearlab.dat
yearlab = yearlab(2:end);
%%
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','--|-|--|:')
  
  
subplot(3,1,1)
plot(yearlab,sigmean(:,1),'LineWidth',0.6);
%hold on;
%plot(yearlab,zeros(length(yearlab)), '-r');
set(gca,'Fontsize',8);
title('GDP per capita (annual %)')
xlim([1954 2017])

subplot(3,1,2)
plot(yearlab,sigmean(:,2),'LineWidth',0.6);
set(gca,'Fontsize',8);
title('TFP')
xlim([1954 2017])

subplot(3,1,3)
plot(yearlab,sigmean(:,3),'LineWidth',0.6);
set(gca,'Fontsize',8);
title('BNDES disbursements to GFCF')
xlim([1954 2017])
saveas(gcf,'vol','epsc')

%%
a2 = zeros(t,2);
a3 = zeros(t,2);
a4 = zeros(t,2);

a6 = zeros(t,2);
a7 = zeros(t,2);
a8 = zeros(t,2);

a10 = zeros(t,2);
a11 = zeros(t,2);
a12 = zeros(t,2);

% Equa��o da produtividade (a2, a3, a4) (produtividade, pib, capital)

for i=1:t
    x=squeeze(Bt_save(2,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a2(i,:)=([LB,UB]);
end
size(a2);

for i=1:t
    x=squeeze(Bt_save(3,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a3(i,:)=([LB,UB]);
end
size(a3);


for i=1:t
    x=squeeze(Bt_save(4,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a4(i,:)=([LB,UB]);
end
size(a4);

% Equa��o do PIB (a6, a7, a8) (produtividade, pib, capital)

for i=1:t
    x=squeeze(Bt_save(6,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a6(i,:)=([LB,UB]);
end
size(a6);

for i=1:t
    x=squeeze(Bt_save(7,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a7(i,:)=([LB,UB]);
end
size(a7);


for i=1:t
    x=squeeze(Bt_save(8,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a8(i,:)=([LB,UB]);
end
size(a8);

% Equa��o do capital (a10, a11, a12) (produtividade, pib, capital)


size(a10);
for i=1:t
    x=squeeze(Bt_save(10,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a10(i,:)=([LB,UB]);
end


for i=1:t
    x=squeeze(Bt_save(11,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a11(i,:)=([LB,UB]);
end

for i=1:t
    x=squeeze(Bt_save(12,i,:));
    [LB,UB] = HPD_SIM(x,0.15);
    a12(i,:)=([LB,UB]);
end
%%
p = [0.16 0.84];
q = [0.10 0.90];
quantil = quantile(Bt_save,p,3);
quantil90 = quantile(Bt_save,q,3);

error3 = squeeze(quantil(3,:,:));
error4 = squeeze(quantil(4,:,:));
error6 = squeeze(quantil(6,:,:));
error7 = squeeze(quantil(7,:,:));
error8 = squeeze(quantil(8,:,:));
error10 = squeeze(quantil(10,:,:));
error11 = squeeze(quantil(11,:,:));
error12 = squeeze(quantil(12,:,:));

error3n = squeeze(quantil90(3,:,:));
error4n = squeeze(quantil90(4,:,:));
error6n = squeeze(quantil90(6,:,:));
error7n = squeeze(quantil90(7,:,:));
error8n = squeeze(quantil90(8,:,:));
error10n = squeeze(quantil90(10,:,:));
error11n = squeeze(quantil90(11,:,:));
error12n = squeeze(quantil90(12,:,:));

%%
figure;
set(0,'DefaultAxesColorOrder','remove') 
plot(yearlab,Bt_postmean(2:4,:));
hold on;
plot(yearlab,Bt_postmean(6:8,:));
hold on;
plot(yearlab,Bt_postmean(10:12,:));
legend('2','3','4','6','7','8','10','11','12')
%%
% TFP no GDP 
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':')
figure;
plot(yearlab,Bt_postmean(3,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error3(:,1)', error3(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error3n(:,1)', error3n(:,2)'  , [0.7 0.8 0.8]);
hold off
alpha(.3)
xlim([1954 2017])
set(gca,'FontSize',12)
%set(gca,'Fontsize',8);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'GDPTFP_1','epsc')
%%
% BNDES no GDP
figure;
plot(yearlab,Bt_postmean(4,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error4(:,1)', error4(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error4n(:,1)', error4n(:,2)'  , [0.7 0.8 0.8]);
hold off
alpha(.3)
%set(gca,'Fontsize',8);
xlim([1954 2017])
set(gca,'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
%title('Efeito das Importa��es na Produtividade ')
saveas(gcf,'GDPBNDES_1','epsc')
%%
% BNDES no TFP
figure;
plot(yearlab,Bt_postmean(8,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error8(:,1)', error8(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error8n(:,1)', error8n(:,2)'  , [0.7 0.8 0.8]);
hold off
alpha(.3)
xlim([1954 2017])
set(gca,'FontSize',12)
%set(gca,'Fontsize',8);
%title('Efeito da Produtividade no Investimento')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'TFPBNDES_1','epsc')

figure;
plot(yearlab,Bt_postmean(12,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error12(:,1)', error12(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error12n(:,1)', error12n(:,2)'  , [0.7 0.8 0.8]);
hold off
alpha(.3)
xlim([1954 2017])
set(gca,'FontSize',12);
%title('Subsidy on Subsidy')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'BNDESBNDES_1','epsc')

figure;
plot(yearlab,Bt_postmean(6,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error6(:,1)', error6(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error6n(:,1)', error6n(:,2)'  , [0.7 0.8 0.8]);
hold off
alpha(.3)
set(gca,'FontSize',12);
%title('Growth on PTF')
xlim([1954 2017])
%set(gca,'FontSize',12)
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'TFPGDP_1','epsc')


%%
if sensitivity==1 
    load('Model2')
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':|-|--|:')
 plot(1:nhor,squeeze(imp75XY(2,1,:)))
    hold on
    plot(1:nhor,squeeze(imp81XY(2,1,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,1,:)))
    legend('1963','1979','2015')
set(gca,'Fontsize',8);
 xlim([1 nhor]);
%title(['Response Growth'],'FontSize',8,'FontWeight','bold');
saveas(gcf,'SResponseGrowth','epsc')


figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':|-|--|:')
 plot(1:nhor,squeeze(imp75XY(2,2,:)))
    hold on
    plot(1:nhor,squeeze(imp81XY(2,2,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,2,:)))
 xlim([1 nhor]);
%title(['Response TFP'],'FontSize',8,'FontWeight','bold');
saveas(gcf,'SResponseTFP','epsc')

figure;
  plot(1:nhor,squeeze(imp75XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(imp81XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,3,:)))
 xlim([1 nhor]);
%title(['Response BNDES'],'FontSize',10,'FontWeight','bold');
saveas(gcf,'SResponseBNDES','epsc')

%%
% TFP no GDP 
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':')
figure;
load('Model1')
plot(yearlab,Bt_postmean(3,:),'-','LineWidth',1);
hold on;
load('Model2')
plot(yearlab,Bt_postmean(3,:),':','LineWidth',1);
alpha(.3)
set(gca,'Fontsize',8);
legend('benchmark prior', 'alternative prior')
xlim([1954 2017])
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'SGDPTFP_1','epsc')
%%
% BNDES no GDP
figure;
load('Model1')
plot(yearlab,Bt_postmean(4,:),'-','LineWidth',1);
hold on;
load('Model2')
plot(yearlab,Bt_postmean(4,:),':','LineWidth',1);
alpha(.3)
set(gca,'Fontsize',8);
legend('benchmark prior', 'alternative prior')
xlim([1954 2017])
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'SGDPBNDES_1','epsc')

figure;
load('Model1')
plot(yearlab,Bt_postmean(8,:),'-','LineWidth',1);
hold on;
load('Model2')
plot(yearlab,Bt_postmean(8,:),':','LineWidth',1);
alpha(.3)
set(gca,'Fontsize',8);
legend('Model1', 'Model2')
xlim([1954 2017])
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'STFPBNDES_1','epsc')





if istore == 1
    load('Model1')

    M1_imp75XY=squeeze(quantile(imp75,qus));
    M1_imp81XY=squeeze(quantile(imp81,qus));
    M1_imp96XY=squeeze(quantile(imp96,qus));
    
     load('Model2')
    M2_imp75XY=squeeze(quantile(imp75,qus));
    M2_imp81XY=squeeze(quantile(imp81,qus));
    M2_imp96XY=squeeze(quantile(imp96,qus));
    %%
figure;

set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|:|--|:')
  subplot(1,3,1);
 plot(1:nhor,squeeze(M1_imp75XY(2,1,:)))
 hold on; 
  plot(1:nhor,squeeze(M2_imp75XY(2,1,:)))

subplot(1,3,2);
    plot(1:nhor,squeeze(M1_imp81XY(2,1,:)))
    hold on;
 plot(1:nhor,squeeze(M2_imp81XY(2,1,:)))
 
   subplot(1,3,3);
    plot(1:nhor,squeeze(M1_imp96XY(2,1,:)))
    hold on;
     plot(1:nhor,squeeze(M2_imp96XY(2,1,:)))
    %legend('1963','1979','2015');
    
    set(gca,'Fontsize',8);
    legend('benchmark prior', 'alternative prior')
 xlim([1 nhor]);
  %ylim([-0.1 0.6]);

%title(['Response Growth'],'FontSize',8,'FontWeight','bold');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
saveas(gcf,'sensiResponseGrowth','epsc')
%%

figure;

set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|:|--|:')
  subplot(1,3,1);
 plot(1:nhor,squeeze(M1_imp75XY(2,2,:)))
 hold on; 
  plot(1:nhor,squeeze(M2_imp75XY(2,2,:)))

subplot(1,3,2);
    plot(1:nhor,squeeze(M1_imp81XY(2,2,:)))
    hold on;
 plot(1:nhor,squeeze(M2_imp81XY(2,2,:)))
 
   subplot(1,3,3);
    plot(1:nhor,squeeze(M1_imp96XY(2,2,:)))
    hold on;
     plot(1:nhor,squeeze(M2_imp96XY(2,2,:)))
    %legend('1963','1979','2015');
    
    set(gca,'Fontsize',8);
    legend('benchmark prior', 'alternative prior')
 xlim([1 nhor]);
  %ylim([-0.002 0.014])%title(['Response TFP'],'FontSize',8,'FontWeight','bold');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 9]);
  saveas(gcf,'sensiResponseTFP','epsc')
  
  
  figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder',':|-|--|:')
 plot(1:nhor,squeeze(M2_imp75XY(2,2,:)))
    hold on
    plot(1:nhor,squeeze(M2_imp81XY(2,2,:)))
    hold on
    plot(1:nhor,squeeze(M2_imp96XY(2,2,:)))
 xlim([1 nhor]);
%title(['Response TFP'],'FontSize',8,'FontWeight','bold');

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 16 6]);
  saveas(gcf,'sensiResponseTFPall','epsc')
 

figure;
  plot(1:nhor,squeeze(M2_imp75XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(M2_imp81XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(M2_imp96XY(2,3,:)))
 xlim([1 nhor]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 16 6]);
saveas(gcf,'sensiResponseBNDESall','epsc')

    
end


end