% plot

figure;
plot(s_grid, pistar_s)

%%
figure;
%subplot(1,2,1)
plot(s_grid,nanmean(nanmean(SUBnbarValues(:,:,:),3)));
xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,nanmean(nanmean(NONnbarValues(:,:,:),3)),'-r');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('earmarked','non-earmarked', 'Location', 'northwest')

figure;
%subplot(1,2,1)
plot(s_grid,nanmean(nanmean(POORnbarValues(:,:,:),3)));
xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,nanmean(nanmean(GOODnbarValues(:,:,:),3)),'-r');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('poor credit','good credit', 'Location', 'northwest')

%%
a11 = squeeze(Params.ebar(1,:,1,1));
a21 = squeeze(Params.ebar(1,:,2,1));
a12 = squeeze(Params.ebar(1,:,1,2));
a22 = squeeze(Params.ebar(1,:,2,2));

figure; 
subplot(2,2,1)
stem(s_grid,a11, 'b')
title('Tax and Good Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Entry Decision') 
subplot(2,2,2)
stem(s_grid,a21, 'b')
title('Sub and Good Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
subplot(2,2,3)
stem(s_grid,a12,'r')
title('Tax and Bad Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Entry Decision') 
subplot(2,2,4)
stem(s_grid,a22, 'r')
title('Sub and Bad Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 

%%

%%
a11 = squeeze(ExitPolicy(:,:,1,1));
a21 = squeeze(ExitPolicy(:,:,2,1));
a12 = squeeze(ExitPolicy(:,:,1,2));
a22 = squeeze(ExitPolicy(:,:,2,2));

figure; 
subplot(2,2,1)
stem3(a11(end-30:end,:), 'Color', 'b','Marker' , 'x')
title('Tax and Good Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
subplot(2,2,2)
stem3(a21(end-30:end,:), 'Color', 'b','Marker' , 'x')
title('Sub and Good Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
zlim([0 1])
subplot(2,2,3)
stem3(a12(end-30:end,:), 'Color', 'r','Marker' , 'x', 'LineStyle', ':')
title('Tax and Bad Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
subplot(2,2,4)
stem3(a22(end-30:end,:), 'Color', 'r','Marker' , 'x', 'LineStyle', ':')
title('Sub and Bad Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 

%%
fprintf('Percentage tax and good credit %8.2f\n', sum(sum(StationaryDist.pdf(:,:,1,1))))

fprintf('Percentage tax and bad credit  %8.2f\n', sum(sum(StationaryDist.pdf(:,:,1,2))))

fprintf('Percentage sub and good credit %8.2f\n', sum(sum(StationaryDist.pdf(:,:,2,1))))

fprintf('Percentage sub and bad credit  %8.2f\n', sum(sum(StationaryDist.pdf(:,:,2,2))))



