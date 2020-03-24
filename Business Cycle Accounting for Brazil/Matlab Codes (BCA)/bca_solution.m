function [gx,hx]=bca_solution(x1)

% This program derives the solution for the Small Open Economy Model

rhoA = x1(1); % AR(1) Productivity
rhoL = x1(2); % AR(1) Labor Wedge
rhoK = x1(3); % AR(1) Capital Wedge
rhoB = x1(4); % AR(1) Bond Wedge

bca_model_run;
