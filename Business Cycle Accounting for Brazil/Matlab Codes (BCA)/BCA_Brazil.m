% This program performs a Business Cycle Accounting excercise for argentina
clear all;
close all;
% Load Data
tic

norm =9 ;

load data_brazil.m;

yt = data_brazil(:,4);
ct = data_brazil(:,1);
ht = data_brazil(:,5);
it = data_brazil(:,2);
kt = data_brazil(:,6);

% Initial Values for Parameters

rhoA = 0.99 ; % AR(1) Productivity
rhoL = 0.99; % AR(1) Labor Wedge
rhoK = 0.99; % AR(1) Capital Wedge
rhoB = 0.99; % AR(1) Bond Wedge
sigA = 0.03; % S.D. Productivity
sigL = 0.03; % S.D. Labor Wedge
sigK = 0.03; % S.D. Capital Wedge
sigB = 0.03; % S.D. Bond Wedge

x0=[rhoA;rhoL;rhoK;rhoB;sigA;sigL;sigK;sigB]; % Initial Values

% 1) Estimate the shocks processes and standard errors

[x1,x1se]=bca(x0,yt,ct,ht,it);

% 2) Solution of the Model

[gx,hx]=bca_solution(x1);

% 3) Compute the wedges 

[at,tht,tkt,tbt]=bca_wedge_alternative2(x1,gx,hx,yt,ct,ht,it,kt);

% 4) Evaluate the shocks into the model


new_graph
new_graph_2
new_graph_3
graph_bca;


% 5) Business Cycle Statistics

%statistics_bca;
toc