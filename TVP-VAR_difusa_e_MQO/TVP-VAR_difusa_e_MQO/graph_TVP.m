%close all;
% Gráficos do Modelo TVP-VAR 
load yearlab.dat
%% Figura das Séries Originais
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
  
subplot(3,1,1)
plot(yearlab,Y(:,1),'-','LineWidth',0.7);
set(gca,'Fontsize',8);
title('BNDES/FBCF')
xlim([1954 2017])

subplot(3,1,2)
plot(yearlab,Y(:,2),'-','LineWidth',0.7);
set(gca,'Fontsize',8);
title('TFP')
xlim([1954 2017])

subplot(3,1,3)
plot(yearlab,Y(:,3),'-','LineWidth',0.7);
set(gca,'Fontsize',8);
title('GDP per capita')
xlim([1954 2017])




%% Desvio Padrão dos resíduos
yearlab = yearlab(p+1:end);
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
  
  
subplot(3,1,1)
plot(yearlab,sigmean(:,1),':','LineWidth',1);
set(gca,'Fontsize',8);
title('BNDES/FBCF')
xlim([1954 2017])

subplot(3,1,2)
plot(yearlab,sigmean(:,2),':','LineWidth',1);
set(gca,'Fontsize',8);
title('TFP')
xlim([1954 2017])

subplot(3,1,3)
plot(yearlab,sigmean(:,3),':','LineWidth',1);
set(gca,'Fontsize',8);
title('GDP per capita')
xlim([1954 2017])


%%
figure;
plot(yearlab,sigmean(:,1),':','LineWidth',1);
set(gca,'Fontsize',8);
hold on
plot(yearlab,sigmean(:,2),'-','LineWidth',1);
hold on;
plot(yearlab,sigmean(:,3),'--','LineWidth',1);
xlim([1954 2017])
legend('BNDES/FBCF','TFP','GDP per capita')

%% Cálculo das Funções Impulso Resposta
% Os anos das Funções estão no script FIR

if istore == 1 
    qus = [.10, .5, .90];
    imp75XY=squeeze(quantile(imp75,qus));
    imp81XY=squeeze(quantile(imp81,qus));
    imp96XY=squeeze(quantile(imp96,qus));
    

 %% Gráfico das Funções Impulso Resposta 
 

 figure;
 set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|:|--')
  
 subplot(3,1,1)
 plot(1:nhor,squeeze(imp75XY(2,1,:)))
     title('Impulso resposta da BNDES/FBCF')
    hold on
    plot(1:nhor,squeeze(imp81XY(2,1,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,1,:)))
    legend('2003','2009','2015');
set(gca,'Fontsize',8);
 xlim([1 nhor]);
title(['Resposta da BNDES/FBCF'],'FontSize',8,'FontWeight','bold');
    
 subplot(3,1,2)
  plot(1:nhor,squeeze(imp75XY(2,2,:)))
    hold on
    plot(1:nhor,squeeze(imp81XY(2,2,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,2,:)))
 xlim([1 nhor]);
title(['Resposta ds TFP'],'FontSize',8,'FontWeight','bold');

subplot(3,1,3)
  plot(1:nhor,squeeze(imp75XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(imp81XY(2,3,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,3,:)))
 xlim([1 nhor]);
title(['Resposta ds TFP'],'FontSize',8,'FontWeight','bold');
       
%%
 % Resposta da Última Equação
    figure;
     set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|:|--')
     plot(1:nhor,squeeze(imp75XY(2,1,:)))
     title('Impulso resposta do GDP per capita')
    hold on
    plot(1:nhor,squeeze(imp81XY(2,1,:)))
    hold on
    plot(1:nhor,squeeze(imp96XY(2,1,:)))
    xlim([1 nhor])
    legend('2003','2009','2015')


%% Intervalos de Confiança para o Parâmetros

a2 = zeros(t,2);
a3 = zeros(t,2);
a4 = zeros(t,2);

a6 = zeros(t,2);
a7 = zeros(t,2);
a8 = zeros(t,2);

a10 = zeros(t,2);
a11 = zeros(t,2);
a12 = zeros(t,2);

% Equação da produtividade (a2, a3, a4) (produtividade, pib, capital)

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

% Equação do PIB (a6, a7, a8) (produtividade, pib, capital)

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

% Equação do capital (a10, a11, a12) (produtividade, pib, capital)


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



           
     %% Figura dos parâmetros variando ao longo do tempo
     
     % Primeira Equação
          figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
subplot(3,1,1)
plot(yearlab,Bt_postmean(2,:),'-','LineWidth',1);
hold on
plot(yearlab,a2,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito da Produtividade na Produtividade ')
xlim([1954 2017])

subplot(3,1,2)
plot(yearlab,Bt_postmean(3,:),'-','LineWidth',1);
hold on
plot(yearlab,a3,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito da Importação na Produtividade')
xlim([1954 2017])

subplot(3,1,3)
plot(yearlab,Bt_postmean(4,:),'-','LineWidth',1);
hold on
plot(yearlab,a4,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito do Investimento na Produtividade')
xlim([1954 2017])

% Segunda Equação 
     figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
subplot(3,1,1)
plot(yearlab,Bt_postmean(6,:),'-','LineWidth',1);
hold on
plot(yearlab,a6,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito da Produtividade nas Importações ')
xlim([1954 2017])

subplot(3,1,2)
plot(yearlab,Bt_postmean(7,:),'-','LineWidth',1);
hold on
plot(yearlab,a7,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito da Importação nas Importações')
xlim([1954 2017])

subplot(3,1,3)
plot(yearlab,Bt_postmean(8,:),'-','LineWidth',1);
hold on
plot(yearlab,a8,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito do Investimento nas Importações')
xlim([1954 2017])


  
  % Terceira Equação 
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
subplot(3,1,1)
plot(yearlab,Bt_postmean(10,:),'-','LineWidth',1);
hold on
plot(yearlab,a10,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito da Produtividade no Investimento')
xlim([1954 2017])

subplot(3,1,2)
plot(yearlab,Bt_postmean(11,:),'-','LineWidth',1);
hold on
plot(yearlab,a11,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito da Importação no Investimento')
xlim([1954 2017])

subplot(3,1,3)
plot(yearlab,Bt_postmean(12,:),'-','LineWidth',1);
hold on
plot(yearlab,a12,':','LineWidth',1);
set(gca,'Fontsize',8);
title('Efeito do Investimento no Investimento')
xlim([1954 2017])



%% Comportamento dos Coeficientes (com quantiles)


p = [0.16 0.84];
q = [0.10 0.90];
quantil = quantile(Bt_save,p,3);
quantil90 = quantile(Bt_save,q,3);

error3 = squeeze(quantil(3,:,:));
error4 = squeeze(quantil(4,:,:));
error6 = squeeze(quantil(6,:,:));
error8 = squeeze(quantil(8,:,:));
error10 = squeeze(quantil(10,:,:));
error12 = squeeze(quantil(12,:,:));

error3n = squeeze(quantil90(3,:,:));
error4n = squeeze(quantil90(4,:,:));
error6n = squeeze(quantil90(6,:,:));
error8n = squeeze(quantil90(8,:,:));
error10n = squeeze(quantil90(10,:,:));
error12n = squeeze(quantil90(12,:,:));

  
% Primeira Equação -- Terceiro Coficiente  
subplot(2,1,1)
plot(yearlab,Bt_postmean(3,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error3(:,1)', error3(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error3n(:,1)', error3n(:,2)'  , [0.7 0.7 1]);
hold off
alpha(.3)
set(gca,'Fontsize',8);
title('TFP on Growth Equation ')
xlim([1954 2017])

% Primeira Equação -- Quarto Coficiente  
subplot(2,1,2)
plot(yearlab,Bt_postmean(4,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error4(:,1)', error4(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error4n(:,1)', error4n(:,2)'  , [0.7 0.7 1]);
hold off
alpha(.3)
set(gca,'Fontsize',8);
title('Subsidy on Growth')
xlim([1954 2017])

 figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
         'DefaultAxesLineStyleOrder','-|-.|-.|:|:')
  
% Segunda Equação -- Quarto Coficiente  
subplot(2,1,1)
plot(yearlab,Bt_postmean(8,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error8(:,1)', error8(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error8n(:,1)', error8n(:,2)'  , [0.7 0.7 1]);
hold off
alpha(.3)
set(gca,'Fontsize',8);
title('Subsidy on TFP ')
xlim([1954 2017])

% Terceira Equação -- Terceiro Coficiente  
subplot(2,1,2)
plot(yearlab,Bt_postmean(10,:),'-','LineWidth',1);
hold on;
shadedplot(yearlab, error10(:,1)', error10(:,2)');
hold off
alpha(.3)
hold on;
shadedplot(yearlab, error10n(:,1)', error10n(:,2)'  , [0.7 0.7 1]);
hold off
alpha(.3)
set(gca,'Fontsize',8);
title('Growth on Subsidy')
xlim([1954 2017])


end
