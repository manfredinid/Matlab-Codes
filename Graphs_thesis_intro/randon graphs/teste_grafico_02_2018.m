% Graficos do ensaio 3
clear all;

load ibc_2018.dat;
load ipea_2018.dat; 

%cria vetor temporal

t1 = datetime(2002,12,1);
meses = t1 + calmonths(1:192);

figure;
plot(meses, ipea_2018/ipea_2018(1))
hold on;
plot(meses, ibc_2018/ibc_2018(1));
legend('Ipea-FBCF','IBC-Br', 'Location', 'Northeast')
