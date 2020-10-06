% plot

figure;
plot(s_grid, pistar_s)

%%
nbarValues=shiftdim(ValuesOnGrid(3,:,:,:,:),1);
normalize_employment=nanmin(nonzeros(nbarValues)); % Normalize so that smallest occouring value of nbar in the baseline is equal to 1.
nbarValues=(nbarValues./(normalize_employment)).*(StationaryDist.pdf~=0);


figure;

subplot(1,2,1)
plot(s_grid,mean(nbarValues(:,:,2,1)),'-b');
xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,mean(nbarValues(:,:,1,2)),'-r');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('earmarked and good credit','non-earmarked and poor credit', 'Location', 'northwest')


subplot(1,2,2)
plot(s_grid,mean(nbarValues(:,:,1,1)),'-b');
xlim([0.9 2.5])
%title('non-earmarked')
%xlabel('productivity')
%ylabel('employees')
%subplot(1,2,2)
hold on;
plot(s_grid,mean(nbarValues(:,:,2,2)),'-r');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
legend('non-earmarked and good credit','earmarked and poor credit', 'Location', 'northwest')
%%

figure;
subplot(1,2,1)
plot(s_grid,mean(nbarValues(:,:,1,1)),'-b');
title('non-earmarked and good credit')
xlabel('productivity')
ylabel('employees')
%ylabel('employees')
%subplot(1,2,2)
subplot(1,2,2)
plot(s_grid,mean(nbarValues(:,:,2,2)),'-r');
%xlim([0.9 2.5])
%title('earmarked')
xlabel('productivity')
ylabel('employees')
title('earmarked and poor credit')


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
stem3(a11(1:50,:), 'Color', 'b','Marker' , 'x')
title('Tax and Good Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
subplot(2,2,2)
stem3(a21(1:50,:), 'Color', 'b','Marker' , 'x')
title('Sub and Good Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
zlim([0 1])
subplot(2,2,3)
stem3(a12(1:50,:), 'Color', 'r','Marker' , 'x', 'LineStyle', ':')
title('Tax and Bad Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 
subplot(2,2,4)
stem3(a22(1:50,:), 'Color', 'r','Marker' , 'x', 'LineStyle', ':')
title('Sub and Bad Credit')
ylabel('Capital') 
xlabel('Productivity')
zlabel('Exit Decision') 

%%
fprintf('Percentage tax and good credit %8.2f\n', sum(sum(StationaryDist.pdf(:,:,1,1))))

fprintf('Percentage tax and bad credit  %8.2f\n', sum(sum(StationaryDist.pdf(:,:,1,2))))

fprintf('Percentage sub and good credit %8.2f\n', sum(sum(StationaryDist.pdf(:,:,2,1))))

fprintf('Percentage sub and bad credit  %8.2f\n', sum(sum(StationaryDist.pdf(:,:,2,2))))



