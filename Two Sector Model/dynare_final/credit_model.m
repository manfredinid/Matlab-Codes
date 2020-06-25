<<<<<<< .mine
%%
%clear all;
close all;

betaa =0.9798;  
varpi = 0.5;
al = 1;
ah =2.75;
sigmaa = 1;
deltaa = 0.0125; 
alphaa = 0.399;

% periods in the plot
n=100;

save parameterfile  betaa varpi al ah sigmaa deltaa alphaa


%pi_p0=1.76;
%pi_g0=1.47;

%pi_pF_obser=1.4;
%pi_gF_obser=0.7;

%pi_pF_simul=0.8;
%pi_gF_simul=0.8;

%%
%pi_p0=1.027/1.008; %2011 2.7 p.m.
%pi_g0=1; % 2011 0.8 p.m.

%pi_pF_obser=1.0356/1.008; %2013 2.36 p.m. % 2017 3.56
%pi_gF_obser= 1; %2013 0.6 p.m.  % 2017 0.79
 
%pi_pF_simul=1;
%pi_gF_simul=1;


pi_p0=(1+0.2142)^(1/4); %2011 2.7 p.m.
pi_g0=(1+0.12)^(1/4); % 2011 0.8 p.m.

pi_pF_obser=(1+0.2216)^(1/4); %2013 2.36 p.m. % 2017 3.56
pi_gF_obser= (1+0.11)^(1/4); %2013 0.6 p.m.  % 2017 0.79
 
pi_pF_simul=(1+0.15)^(1/4);%(1.135+1.2657)/2;
pi_gF_simul=(1+0.15)^(1/4);%(1.135+1.2657)/2;

%%
dynare model_observed
load('model_observed_results.mat')
y_obser = y(1:n);
i_obser = i(1:n);
k_obser = k(1:n);
c_obser = c(1:n);
kl_obser = kl(1:n);
kh_obser = kh(1:n);

%% 
dynare model_simul
load('model_simul_results.mat')
y_simul = y(1:n);
i_simul = i(1:n);
k_simul = k(1:n);
c_simul = c(1:n);
kl_simul = kl(1:n);
kh_simul = kh(1:n);


%% Graphs

figure(1)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',0.8);
subplot(2,2,1);
plot(c_obser);
hold on
plot(c_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,2);
plot(i_obser);
hold on
plot(i_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,3);
plot(k_obser);
hold on
plot(k_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,4);
plot(y_obser);
hold on
plot(y_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')

snapnow

% capital-output ratio
figure(2)
plot(k_obser./y_obser);
hold on
plot(k_simul./y_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (observed)'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,1,1);
plot(kl_obser);
hold on
plot(kl_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['low-tecnhology capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')
subplot(2,1,2);
plot(kh_obser);
hold on
plot(kh_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['high-technology capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')
snapnow

%% Selected Graphs

figure(3);
plot(kl_obser/kl_obser(1)*100);
hold on
plot(kl_simul/kl_simul(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['low-tecnhology capital'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
saveas(gcf,'kl.eps','epsc')

figure(4);
plot(kh_obser/kh_obser(1)*100);
hold on
plot(kh_simul/kh_simul(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['high-technology capital'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
saveas(gcf,'kh.eps','epsc')

figure(5);
plot(i_obser/i_obser(1)*100);
hold on
plot(i_simul/i_simul(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
legend({'observed','alternative', 'Location', 'Best'},'FontSize',12)
saveas(gcf,'investment.eps','epsc')

figure(6);
plot(y_obser/y_obser(1)*100);
hold on
plot(y_simul/y_obser(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
saveas(gcf,'output.eps','epsc')
=======
%%
%clear all;
close all;

betaa =0.9;  
varpi = 0.5;
al = 1;
ah =2.75;
sigmaa = 0.63;
deltaa = 0.04; 
alphaa = 0.6;

% periods in the plot
n=40;

save parameterfile  betaa varpi al ah sigmaa deltaa alphaa


%pi_p0=1.76;
%pi_g0=1.47;

%pi_pF_obser=1.4;
%pi_gF_obser=0.7;

%pi_pF_simul=0.8;
%pi_gF_simul=0.8;

%%
%pi_p0=1.027/1.008; %2011 2.7 p.m.
%pi_g0=1; % 2011 0.8 p.m.

%pi_pF_obser=1.0356/1.008; %2013 2.36 p.m. % 2017 3.56
%pi_gF_obser= 1; %2013 0.6 p.m.  % 2017 0.79
 
%pi_pF_simul=1;
%pi_gF_simul=1;


pi_p0=(1+0.2455)^(1/4); %2011 2.7 p.m.
pi_g0=(1+0.1225)^(1/4); % 2011 0.8 p.m.

pi_pF_obser=(1+0.2657)^(1/4); %2013 2.36 p.m. % 2017 3.56
pi_gF_obser= (1+0.1350)^(1/4); %2013 0.6 p.m.  % 2017 0.79
 
pi_pF_simul=(1+0.184)^(1/4);%(1.135+1.2657)/2;
pi_gF_simul=(1+0.184)^(1/4);%(1.135+1.2657)/2;

%%
dynare model_observed
load('model_observed_results.mat')
y_obser = y(1:n);
i_obser = i(1:n);
k_obser = k(1:n);
c_obser = c(1:n);
kl_obser = kl(1:n);
kh_obser = kh(1:n);

%% 
dynare model_simul
load('model_simul_results.mat')
y_simul = y(1:n);
i_simul = i(1:n);
k_simul = k(1:n);
c_simul = c(1:n);
kl_simul = kl(1:n);
kh_simul = kh(1:n);


%% Graphs

figure(1)
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:','DefaultLineLineWidth',0.8);
subplot(2,2,1);
plot(c_obser);
hold on
plot(c_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Consumption'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,2);
plot(i_obser);
hold on
plot(i_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,3);
plot(k_obser);
hold on
plot(k_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,2,4);
plot(y_obser);
hold on
plot(y_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')

snapnow

% capital-output ratio
figure(2)
plot(k_obser./y_obser);
hold on
plot(k_simul./y_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Capital Output Ratio (observed)'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')


subplot(2,1,1);
plot(kl_obser);
hold on
plot(kl_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['low-tecnhology capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')
subplot(2,1,2);
plot(kh_obser);
hold on
plot(kh_simul,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['high-technology capital'],'FontSize',8,'FontWeight','bold');
legend('observed','alternative', 'Location', 'Best')
snapnow

%% Selected Graphs

figure(3);
plot(kl_obser/kl_obser(1)*100);
hold on
plot(kl_simul/kl_simul(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['low-tecnhology capital'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
saveas(gcf,'kl.eps','epsc')

figure(4);
plot(kh_obser/kh_obser(1)*100);
hold on
plot(kh_simul/kh_simul(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['high-technology capital'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
saveas(gcf,'kh.eps','epsc')

figure(5);
plot(i_obser/i_obser(1)*100);
hold on
plot(i_simul/i_simul(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Investment'],'FontSize',8,'FontWeight','bold');
legend({'observed','alternative', 'Location', 'Best'},'FontSize',12)
saveas(gcf,'investment.eps','epsc')

figure(6);
plot(y_obser/y_obser(1)*100);
hold on
plot(y_simul/y_obser(1)*100,'r');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
title(['Output'],'FontSize',8,'FontWeight','bold');
%legend('observed','alternative', 'Location', 'Best')
saveas(gcf,'output.eps','epsc')
>>>>>>> .r181
