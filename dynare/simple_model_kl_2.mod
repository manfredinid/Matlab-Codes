/* Simple RBC compute the steady state values for a credit shock
* This model is based on Landon-Lane and Robertson (2005)
* The model has a two firms types, that differ in productivitty
* The exogenous variables are the total credit in the economy and
* the public credit. The private credit is the residual between total
* and public credit.
*/

//% Variable names
var kh kl c y k i;

// list of exogenous variables 
varexo pi_g pi_p;

// list of parameters
parameters varpi $\varpi$, 
           beta $\beta$,
           al $a_l$,

           ah $a_h$,
           sigma $\sigma$,
           delta $\delta$
           a $(a_h^\varpi)(a_l^{1-\varpi})$,
           alpha $\alpha$,
           theta  $\dfrac{(1-\varpi)}{\varpi}$;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
// parameters values
%Parameters RBC Cycles
beta =0.92;  % ok
varpi = 1/2;
al =  5.5;
ah =5.8;
sigma = 2;
delta = 0.03; % Depreciation rate 
a = (ah^varpi)*(al^(1-varpi));
alpha =0.7;
theta = (1-varpi)/varpi;



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;

%y_obs=log( y )−log( y(−1) )

# psi = ((pi_p*(1-varpi))/(pi_g*varpi));

# A = (ah^varpi)*(al^(1-varpi));

# expH = alpha*varpi;

# expL = alpha*(1-varpi);


[name='Aggregate Output']
exp(y) = A*(exp(kh(-1))^expH)*((psi*exp(kh(-1)))^expL);


[name='Euler Equation']

exp(c)^(-sigma) = beta*(exp(c(+1))^(-sigma))*( (expH*(A*(exp(kh)^expH)*((psi*exp(kh))^expL)))/(exp(kh)*(pi_p+pi_g)) + (expL*(A*(exp(kh)^expH)*((psi*exp(kh))^expL)))/((pi_p+pi_g)*(psi*exp(kh))) + (1-delta));


[name='Budget Constrain']

exp(c)= (A*(exp(kh(-1))^expH)*((psi*exp(kh(-1)))^expL)) - pi_g*(psi*exp(kh) - (1-delta)*(psi*exp(kh(-1)))) - pi_p*(exp(kh) - (1-delta)*exp(kh(-1)));


[name='low-tech capital']
exp(kl) = ((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh);

[name='total capital']
exp(k)= exp(kh) + exp(kl);

exp(i)= exp(y)-exp(c);
end;

%----------------------------------------------------------------
% 4. steadexp(y) state
%----------------------------------------------------------------
steady_state_model;

kh = log((a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alpha*(1-varpi)))/(((1/beta-1+delta)*(pi_p))/(alpha*varpi)))^(1/(1-alpha)));


y = log(a*(exp(kh)^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh))^(alpha*(1-varpi)))); // output

c= log((a*(exp(kh)^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh))^(alpha*(1-varpi)))) - pi_g*delta*(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh)) - pi_p*delta*exp(kh)); // consumption

kl = log(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh));


k= log(exp(kh) + exp(kl));

i= log(exp(y)-exp(c));

end;



%----------------------------------------------------------------
% 5. initial value
%----------------------------------------------------------------
initval;

pi_p=3.51;
pi_g=1;

kl = log(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh));
k= log(exp(kh) + exp(kl));
kh = log((a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alpha*(1-varpi)))/(((1/beta-1+delta)*(pi_p))/(alpha*varpi)))^(1/(1-alpha)));

y = 9.87; // output 2011

c = 9.64; // consumption 2011

i = 8.29; // investimento 2011

end;
//steady(solve_algo=4,maxit=100);
check;

%----------------------------------------------------------------
% 6. end value (permanet shock)
%----------------------------------------------------------------
endval;


pi_g = 1;
pi_p =4.2;


kh = log((a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alpha*(1-varpi)))/(((1/beta-1+delta)*(pi_p))/(alpha*varpi)))^(1/(1-alpha)));


y = log(a*(exp(kh)^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh))^(alpha*(1-varpi)))); // output

c= log((a*(exp(kh)^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh))^(alpha*(1-varpi)))) - pi_g*delta*(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh)) - pi_p*delta*exp(kh)); // consumption

kl = log(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh));


k= 10.9;

i= log(exp(y)-exp(c));

end;
//steady(solve_algo=4,maxit=300);

//figure; plot(log(linspace(exp(3.5),exp(4.12),100)))
//shocks;

//var pi_p;
//periods 1:50;
//values (log(linspace(exp(1.12),exp(1.09),50)));


//var pi_p;
//periods 1:180;
//values (log(linspace(exp(3.6),exp(4.2),180)));
//end;

%----------------------------------------------------------------
% 7. trajectory simulation
%----------------------------------------------------------------
// numerical simulation of the trajectory
// must specify a large number of periods 
simul(periods=180);
perfect_foresight_setup(periods=180);
perfect_foresight_solver(maxit=20,no_homotopy);
//,no_homotopy

resid;

///model_diagnostics(M_, options_, oo_)
%----------------------------------------------------------------
% 8. plots
%----------------------------------------------------------------

//rplot kh;
//rplot kl;
//rplot k;
//rplot y;
//rplot c;
//rplot pi_g;
//rplot pi_p;
//rplot i;

//write_latex_dynamic_model;
//write_latex_static_model;