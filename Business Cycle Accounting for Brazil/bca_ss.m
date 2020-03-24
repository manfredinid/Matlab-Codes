function [b_y, delta, i_y, L, eta, RBAR, alph, nu, gama, n, A, TAUL, TAUK, TAUB, pssi, betta, phi, B,...
        c, cp, ik, ikp, la, lap, y, yp, k, kp, b, bp, r, rp, a, ap, tau_l, tau_lp, tau_k, tau_kp, tau_b, tau_bp]=bca_ss


% Initial parameters for the calibration

b_y=0.3323; % ratio debt/gdp 1995 - 2018 
delta=0.05; % Depreciation Rate
i_y=0.19; % Investment Share
L=0.21; %Labor Supply
eta=0.0001; % Elasticity of bonds
RBAR=1.02; % Interest Rate ROW.
alph=0.3; % Capital Share in Production Function
nu=1.6; % Labor Supply Elasticity
gama=0.0063; % Technological Progress
n=0.0180; % Population Growth
A=1;
TAUL=1;
TAUK=1;


param0=[b_y;i_y;RBAR;alph;delta;L;nu;gama;n];
opciones = optimset('Display','off','TolFun',1e-30,'MaxIter',10000,'TolX',1e-30,'MaxFunEvals',10000);
xg = [0.1; 0.1; 1;0.99;2]; %initial guess
x01 = fsolve(@ss_bca,xg,opciones,param0);

C=x01(1,1); % Steady State Consumption
K=x01(2,1); % Steady State Capital
TAUB=x01(3,1); % Steady State Bond Wedge
BETTA=x01(4,1); % Discount Factor
PSSI=x01(5,1); % Weight on Leisure

Y = ( K ^ alph ) * ( L ^ ( 1 - alph ) ) ; % Production Function
B = b_y * Y; % Debt
I = i_y * Y; % Investment
phi= 0.25 / ( n + gama + gama * n + delta ); % Adjustment Cost of Capital BGG (1999)

% Parameters

pssi=PSSI;
betta=BETTA;

% Variables

c = log(C); % Consumption
ik = log(I); % Investment
la = log(L); % Labor
y = log(Y); % Output
k = log(K); % Capital
b = log(B); % Debt
r = log(RBAR); % Interest Rate
a = log(A); % Productivity
tau_l = log(TAUL); % Labor Wedge
tau_k = log(TAUK); % Capital Wedge
tau_b = log(TAUB); % Bond Wedge

% Variables one period ahead

cp = log(C);
ikp = log(I);
lap = log(L);
yp = log(Y);
kp = log(K);
bp = log(B);
rp = log(RBAR);
ap = log(A);
tau_lp = log(TAUL);
tau_kp = log(TAUK);
tau_bp = log(TAUB);
