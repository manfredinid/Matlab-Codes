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
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 10]);
saveas(gcf,'credit_real','epsc')

%%
load Private
quarter = datetime(2003,01,30):calquarters(1):datetime(2019,12,31);
dateshift(quarter,'end','month');

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');
  plot(quarter, Privatenonfinancial1);
ylabel('Credit-to-GDP ratio')
legend('Brazil', 'Japan', 'Portugal','Location', 'northeast');

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 10]);
saveas(gcf,'Private_non_finantial','epsc')
%%
%%% Import the data
[~, ~, raw] = xlsread('BNDES.xlsx','Sheet4','A2:D65');

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
YEAR = data(:,1);
PIBPC = data(:,2);
TFP = data(:,3);
BNDESFBCF = data(:,4);

% Clear temporary variables
clearvars data raw;
figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|:|--|:');
plot(YEAR, BNDESFBCF);
xlim([1954 2017])
hold on;
line([1974 1974], get(gca, 'ylim'), 'Color', 'red',...
   'LineStyle', ':');
line([1979 1979], get(gca, 'ylim'), 'Color', 'red',...
   'LineStyle', ':');
line([2009 2009], get(gca, 'ylim'), 'Color', 'red',...
   'LineStyle', ':');
line([2014 2014], get(gca, 'ylim'), 'Color', 'red',...
   'LineStyle', ':');
line([1930 1930], get(gca, 'ylim'), 'Color', 'red',...
   'LineStyle', ':');
line([1964 1964], get(gca, 'ylim'), 'Color', 'red',...
   'LineStyle', ':');
txt = {'Authoritarian national','developmentalism'};
text(1974, 23 ,txt,'FontSize',14)
txt = {'New','developmentalism'};
text(2010, 23 ,txt,'FontSize',14,'HorizontalAlignment','center')
txt = {'Golden age of:','import substitution'};
text(1955, 23 ,txt,'FontSize',14)
ylabel('BNDES disbursements to GFCF')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 25 15]);
saveas(gcf,'BNDES','epsc')

%%


set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');
figure;
indicator = [ones(1,sum(YEAR<1964)) zeros(1,sum((YEAR>=1964).*(YEAR<1974)))...
    ones(1,sum((YEAR>=1974).*(YEAR<1980)))  zeros(1,sum((YEAR>=1980).*(YEAR<2009)))...
    ones(1,sum((YEAR>=2009).*(YEAR<2015)))  zeros(1,sum((YEAR>=2015).*(YEAR<=2017)))]';
shadedTimeSeries(YEAR, BNDESFBCF ,[indicator], 'Ano', {'invest'}, [0 0 0]+0.8, 07);
txt = {'\leftarrow Authoritarian national','developmentalism'};
text(1976, 15 ,txt,'FontSize',12)
txt = {'New','developmentalism \rightarrow '};
text(2003, 1,txt,'FontSize',12,'HorizontalAlignment', 'center')
txt = {'\leftarrow Golden age of','import substitution'};
text(1960, 20 ,txt,'FontSize',12)
ylabel('BNDES disbursements to GFCF')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 15]);
saveas(gcf,'BNDES1','epsc')
%%
format shortG
load CF
year = 1995:2019;


figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','--|-.|--|:');
  plot(year, CF);
  hold on
  plot(year, zeros(1,25), '-r')
  xlim([1995 2019])
ylabel('Financial Flows (US$ mill.)')
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 10]);
saveas(gcf,'ck','epsc')

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','--|-.|--|:');
  plot(year(10:end-2), CF(10:end-2));
  hold on
  plot(year(10:end-2), zeros(1,14), '-r')
  xlim([2004 2017])
ylabel('Financial Flows (US$ mill.)')
set(gca, 'YTickLabel', num2cell(get(gca, 'YTick')))

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 10]);
saveas(gcf,'ck2','epsc')
%%
load('credit_real')
% create quarter from July 1989 until April 2016
quarter = datetime(1952,01,30):calquarters(1):datetime(2019,12,31);
quarter =dateshift(quarter,'end','month');

figure;
set(groot,'DefaultAxesColorOrder',[0 0 0],...
      'DefaultAxesLineStyleOrder','-|-.|--|:');
  plot(quarter, credit_real);
  xtickformat('QQQ-yyyy')
  xtickangle(45)
  hold on;
  line([quarter(1) quarter(end)],[0,0],'Color','red','LineStyle','-');
ylabel('Real Credit Growth (Y/Y)')
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 10]);
saveas(gcf,'credit_real','epsc')
