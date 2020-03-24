% This program graph all the results of the model

% 1) Plot wedges (1990s) t=1 (1991) t=16 (2006) t=7 (1997)

At=exp(at);
At=At / At(7) * 100;

THt=exp(tht);
THt= THt / THt(7) * 100;

TKt=exp(tkt);
TKt=TKt / TKt(7) * 100;

TBt=exp(tbt);
TBt=TBt / TBt(7) * 100;

t=1995:2018;
T=length(yt);

figure(1)
%set(gcf,'Color',[1,1,1]);
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:')
  
subplot(2,2,1);
plot(t, At ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Efficiency Wedge (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Produtividade'],'FontSize',8,'FontWeight','bold');

subplot(2,2,2);
plot(t, THt ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Labor Wedge (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Desvio do Trabalho'],'FontSize',8,'FontWeight','bold');

subplot(2,2,3);
plot(t, TKt ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Capital Wedge (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Desvio do Capital'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, TBt ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Bond Wedge (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Desvio dos Títulos'],'FontSize',8,'FontWeight','bold');

% 2) Plot data (1991s) t=1 (1974) t=33 (2006) t=17 (1991) t=22 (1995)

ytg=yt / yt(7) * 100;
ctg=ct / ct(7) * 100;
htg=ht / ht(7) * 100;
itg=it / it(7) * 100;

figure(2)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ytg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['GDP (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['PIB'],'FontSize',8,'FontWeight','bold');

subplot(2,2,2);
plot(t, ctg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Consumption (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Consumo'],'FontSize',8,'FontWeight','bold');

subplot(2,2,3);
plot(t, htg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Labor (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Trabalho'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, itg ,'-','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Investment (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Investimento'],'FontSize',8,'FontWeight','bold');

% 3) Forecast Predictions of the model (1991s) t=1 (1974) t=33 (2006) t=17 (1991) t=22 (1995)

% 3.1) Productivity

xt = [zeros(3,T); at'; zeros(3,T)];

xt = [xt, zeros(7,1)];

yt_at = [zeros(1,T)];

zt = [zeros(4,T)];

for kk=1:T
    
    zt(:,kk) = gx * xt(:,kk);
    
    xtp = hx * xt(:,kk);
    
    xt(1:3,kk+1) = xtp(1:3);
    
    yt_at(kk) = zt(4,kk);
        
end

yt_at = exp(yt_at);
yt_at_bra = yt_at;
yt_at = yt_at / yt_at(7) * 100;

% 3.2) Labor Wedge

xt = [zeros(3,T); zeros(1,T); tht'; zeros(2,T)];

xt = [xt, zeros(7,1)];

yt_tht = [zeros(1,T)];

zt = [zeros(4,T)];

for kk=1:T
    
    zt(:,kk) = gx * xt(:,kk);
    
    %xtp = hx * xt(:,kk);
    
    %xt(1:3,kk+1) = xtp(1:3);
    
    yt_tht(kk) = zt(4,kk);
        
end

yt_tht = exp(yt_tht);
yt_tht_bra = yt_tht;
yt_tht = yt_tht / yt_tht(7) * 100;


% 3.3) Capital Wedge

xt = [zeros(3,T); zeros(2,T); tkt'; zeros(1,T)];

xt = [xt, zeros(7,1)];

yt_tkt = [zeros(1,T)];

zt = [zeros(4,T)];

for kk=1:T
    
    zt(:,kk) = gx * xt(:,kk);
    
    xtp = hx * xt(:,kk);
    
    xt(1:3,kk+1) = xtp(1:3);
    
    yt_tkt(kk) = zt(4,kk);
        
end

yt_tkt = exp(yt_tkt);
yt_tkt_bra = yt_tkt;
yt_tkt = yt_tkt / yt_tkt(7) * 100;


% 3.4) Bond Wedge

xt = [zeros(3,T); zeros(3,T); tbt'];

xt = [xt, zeros(7,1)];

yt_tbt = [zeros(1,T)];

zt = [zeros(4,T)];

for kk=1:T
    
    zt(:,kk) = gx * xt(:,kk);
    
    xtp = hx * xt(:,kk);
    
    xt(1:3,kk+1) = xtp(1:3);
    
    yt_tbt(kk) = zt(4,kk);
        
end

yt_tbt = exp(yt_tbt);
yt_tbt_bra = yt_tbt;
yt_tbt = yt_tbt / yt_tbt(7) * 100;

%---------TESTE
% Aqui a matriz tem que ter dimensão 7. Então é só ligar os wedges que
% deseja e cuidar para a dimensão ser 7 com os zeros da matriz
xt = [zeros(5,T); at'; tht'];

xt = [xt, zeros(7,1)];

yt_that = [zeros(1,T)];

zt = [zeros(4,T)];

for kk=1:T
    
    zt(:,kk) = gx * xt(:,kk);
    
    xtp = hx * xt(:,kk);
    
    xt(1:3,kk+1) = xtp(1:3);
    
    yt_that(kk) = zt(4,kk);
    
    yt_that(kk) = zt(4,kk);
        
end

yt_that = exp(yt_that);
yt_that_bra = yt_that;
yt_that = yt_that / yt_that(7) * 100;
%---------TESTE
%---------TESTE
xt = [zeros(3,T); at'; zeros(3,T)];

xt = [xt, zeros(7,1)];

ht_tht = [zeros(1,T)];

zt = [zeros(4,T)];

for kk=1:T
    
    zt(:,kk) = gx * xt(:,kk);
    
    xtp = hx * xt(:,kk);
    
    xt(1:3,kk+1) = xtp(1:3);
    
    ht_tht(kk) = zt(3,kk);
        
end

ht_tht = exp(ht_tht);
ht_tht_bra = ht_tht;
ht_tht = ht_tht / ht_tht(7) * 100;
%---------TESTE


figure(3)
set(gcf,'Color',[1,1,1]);

subplot(2,2,1);
plot(t, ytg,'-',t, yt_at,':','LineWidth',1);
%legend('Data','Model',2);
legend('Dados','Modelo',2);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Effect of Productivity on GDP (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Efeito da Produtividade no PIB'],'FontSize',8,'FontWeight','bold');

subplot(2,2,2);
plot(t, ytg,'-',t, yt_tht,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Effect of Labor Wedge on GDP (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Efeito do Desvio do Trabalho no PIB'],'FontSize',8,'FontWeight','bold');

subplot(2,2,3);
plot(t, ytg,'-',t, yt_tkt,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Effect of Capital Wedge on GDP (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Efeito do Desvio do Capital no PIB'],'FontSize',8,'FontWeight','bold');

subplot(2,2,4);
plot(t, ytg,'-',t, yt_tbt,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Effect of Bond Wedge on GDP (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Efeito do Desvio dos Títulos no PIB'],'FontSize',8,'FontWeight','bold');


figure;
plot(t, ytg,'-.',t, yt_that,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Effect of TFP and Labor Wedge on GDP (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Teste'],'FontSize',8,'FontWeight','bold');

figure;
plot(t, htg,'-.',t, ht_tht,':','LineWidth',1);
set(gca,'Fontsize',8);
xlim([1995 2018]);
%title(['Effect of Labor Wedge on Labor (Brazil)'],'FontSize',8,'FontWeight','bold');
title(['Teste'],'FontSize',8,'FontWeight','bold');

