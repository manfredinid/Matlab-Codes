%Gráficos da Tese - PPGECO
clear all;
%Abril de 2018

load Gross_capital_formation.dat;
load anual_1960_2018.dat; 
load trimestres.dat;
load investimento_trimestral.dat;
load trimestre.dat;
load producao_industrial_2012_100.dat; % [Bens de capital Bens intermediários Bens de consumo]
load FBKF_IPEA.dat; % [Construção Civil	Consumo Aparente de Máquinas e Equipamentos	Indicador Ipea de FBCF]
load ibcbr.dat;
load pqBI.dat;

anos=anual_1960_2018; %vetor de anos


%meses de janeiro de 2000 até janeiro de 2018
t1 = datetime(1995,12,1);
meses = t1 + calmonths(1:265);

%trimestres 1996-Q1 até 2018-Q1
trimestres=trimestre;



%Formação Bruta de Capital (%PIB) - Dados do Banco Mundial
GKF_GDP=Gross_capital_formation;
ano=anos(1:57);


            
%Indicador anual do Milagre, II PND, Anos 1980 e Espetáculo do crescimento;
%1965-1973 1974-1979 1980-1989 2008-2015

figure;
indicator = [zeros(1,5) ones(1,26) zeros(1,16) ones(1, 8) zeros(1,2)]';
shadedTimeSeries( ano,GKF_GDP ,[indicator], 'Ano', {'invest'}, [0 0 0]+0.8, 07);



%Taxa de investimento a preços correntes, obtida a partir da relação entre a Formação bruta de capital fixo e o Produto interno bruto trimestral nominal (IBGE)
%Dados IBGE de PIB a preços de mercado e Formação Bruta de Capital Fixo
trimestre=trimestres(1:88);
FBKF_Q=investimento_trimestral;

figure;
indicator1 = [zeros(1,44) ones(1,30) zeros(1,14)]';
shadedTimeSeries( trimestre,FBKF_Q ,[indicator1], 'Trimestre', {'invest'}, [0 0 0]+0.8);

%investimento e ciclo econômico

figure;
%findpeaks(FBKF_Q,trimestre); %Aponta os picos da série
%hold on;
indicatorciclo = [zeros(1,8)  ones(1,5) zeros(1,8) ones(1,3) zeros(1,4) ones(1,2) zeros(1,21) ones(1,2) zeros(1,20) ones(1,11) zeros(1,4)]';
shadedTimeSeries( trimestre,FBKF_Q ,[indicatorciclo], 'Trimestre', {'invest'}, [0 0 0]+0.8);




%Indicador IPEA de formação bruta de capital fixo
%Janeiro de 1996 até agosto de 2017
mes=meses(1:260);
FBKF_IPEA_I=FBKF_IPEA(:,3);
figure;
plot(mes, FBKF_IPEA);


%Indicador IPEA com datação do ciclo econômico
mes = datenum(mes);
indicatorciclomes = [zeros(1,24)  ones(1,15) zeros(1,24) ones(1,9) zeros(1,12) ones(1,6) zeros(1,63) ones(1,6) zeros(1,60) ones(1,33) zeros(1,8)]';
figure;
shadedTimeSeries(mes,FBKF_IPEA_I,[indicatorciclomes], 'mês', {'invest'}, [0 0 0]+0.8);
datetick('x','mm-yyyy','keeplimits','keepticks');


mes = datenum(mes);
indicatorciclomes = [zeros(1,24)  ones(1,15) zeros(1,24) ones(1,9) zeros(1,12) ones(1,6) zeros(1,63) ones(1,6) zeros(1,60) ones(1,33) zeros(1,8)]';
figure;
shadedTimeSeries(mes,FBKF_IPEA_I,[indicatorciclomes], 'mês', {'invest'}, [0 0 0]+0.8);
datetick('x','mm-yyyy','keeplimits','keepticks');



%IBC-br Janeiro início em janeiro de 2003 
mes=meses(85:260);
FBKF_IPEA_I=FBKF_IPEA(85:260,3);
ibcbr=ibcbr;
S=[FBKF_IPEA_I ibcbr];


indicatorciclomes = [ones(1,6) zeros(1,63) ones(1,6) zeros(1,60) ones(1,33) zeros(1,8)]';
mes = datenum(mes);
figure;
shadedTimeSeries(mes,FBKF_IPEA_I,[indicatorciclomes], 'mês', {}, [0 0 0]+0.8);
datetick('x','mm-yyyy','keeplimits','keepticks');
hold on;
plot(mes, ibcbr);


figure;
y = unnamed;
x = unnamed2;
yyaxis left
plot(unnamed1, y);
hold on;
yyaxis right
plot(unnamed1, x);


%cria vetor temporal

series= [FBKF_IPEA]
ts1=timeseries(series)
ts1.TimeInfo.StartDate = 'Jan-2002';     % Set start date.
ts1.TimeInfo.Format = 'mm, yyyy';       % Set format for display on x-axis.


load dados.mat;
anos=anual_1960_2018; %vetor de anos
figure;
plot(dados(1,:))
