% Graphs

%% Load data
clear all;
close all;


load ICE
load GAD

% Rename 

Current = ndicedaSituaoAtual;

Expected = ndicedeClimaEconmico;

ICE = ndicedeExpectativas;

% create quarter from July 1989 until April 2016
quarter = datetime(1989,07,30):calquarters(1):datetime(2016,04,31);
dateshift(quarter,'end','month');

% Start in 2001
quarter_2001 = quarter(55:end);
ICE_2001 = ICE(55:end);
Expected_2001 = Expected(55:end);

% Start in 1995
quarter_1995 = quarter(23:end);
ICE_1995 = ICE(23:end);
Expected_1995 = Expected(23:end);

% Another way to create quarters
quarter_gad = datetime(2000,01,30):calquarters(1):datetime(2015,12,31);
dateshift(quarter_gad,'end','month');

%% Economic Climate Index (ECI)
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-.|-|--|:');
plot(quarter, ICE);
%hold on;
%plot(quarter, Expected);
hold on;
line([quarter(1) quarter(end)],[100,100],'Color','red','LineStyle','-');
set(gca,'Fontsize',8);
%xlim([2003 2017]);
hold off;
%legend('ICE','Expected');
snapnow

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');
plot(quarter_2001, ICE_2001);
hold on;
plot(quarter_2001, Expected_2001);
hold on;
ylim([30 180]);
line([quarter_2001(1) quarter_2001(end)],[100,100],'Color','red','LineStyle','-');
legend('Business Confidence Index (ICE)','Business Expectations Index (IE)','Location', 'southwest');
snapnow

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');
plot(quarter_1995, ICE_1995);
hold on;
plot(quarter_1995, Expected_1995);
hold on;
ylim([30 180]);
line([quarter_1995(1) quarter_1995(end)],[100,100],'Color','red','LineStyle','-');
legend('Business Confidence Index (ICE)','Business Expectations Index (IE)','Location', 'southwest');
snapnow

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');
plot(quarter_gad, GAD);
%hold on;
%plot(quarter_1995, Expected_1995);
hold on;
%ylim([30 180]);
line([quarter_gad(1) quarter_gad(end)],[mean(GAD),mean(GAD)],'Color','red','LineStyle','-');
legend('Antidumping Investigations', 'Average','Location', 'northwest');
snapnow