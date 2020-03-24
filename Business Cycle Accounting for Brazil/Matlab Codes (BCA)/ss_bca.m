function q = ss_bca(x,param0)

% param0=[b_y;i_y;RBAR;alph;delta;L;nu;gama;n];

b_y=param0(1,1); i_y=param0(2,1); RBAR=param0(3,1); alph=param0(4,1); delta=param0(5,1);
L=param0(6,1); nu=param0(7,1); gama=param0(8,1); n=param0(9,1);

C=x(1,1); K=x(2,1); TAUB=x(3,1); BETTA=x(4,1); PSSI=x(5,1);

% Derivatives Utility Function

UC = 1 / C; % Marginal Utility of Consumption 
UL = -PSSI * inv( 1 - L ); % Marginal Disutility of Labor

% Derivatives of Production Function

F_L = (K^alph) * (L^(-alph)) * (1-alph) ; %Marginal Productivity of Labor
F_K = (K^(alph-1)) * (L^(1-alph)) * alph; %Marginal Productivity of Capital

% Some Definitions

Y = ( K ^ alph ) * ( L ^ ( 1 - alph ) ) ; % Production Function
B = b_y * Y; % Debt

q(1) = F_L + UL / UC ; % Consumption - Leisure condition
q(2) = 1 - BETTA / (1 + gama) * ( F_K + (1 - delta) ); % Euler Equation Capital
q(3) = 1 - BETTA / (1 + gama) * TAUB * RBAR ; % Euler Equation Bond
q(4) = B * (1 + gama + n + gama * n - RBAR) - C - i_y * Y + Y; % Resource Constraint
q(5) = i_y * Y / K - ( n + gama + gama * n + delta ); % Capital Accumulation