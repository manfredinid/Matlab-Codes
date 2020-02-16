/* Simple RBC compute the steady state values for a credit shock
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
al =  0.5;
ah =0.8;
sigma = 2;
delta = 0.03; % Depreciation rate 
a = (ah^varpi)*(al^(1-varpi));
alpha =0.5;
theta = (1-varpi)/varpi;



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;


# psi = ((pi_p*(1-varpi))/(pi_g*varpi));

# A = (ah^varpi)*(al^(1-varpi));

# expH = alpha*varpi;

# expL = alpha*(1-varpi);


[name='Aggregate Output']
y = A*(kh(-1)^expH)*((psi*kh(-1))^expL);


[name='Euler Equation']

c^(-sigma) = beta*(c(+1)^(-sigma))*( (expH*(A*(kh^expH)*((psi*kh)^expL)))/(kh*(pi_p+pi_g)) + (expL*(A*(kh^expH)*((psi*kh)^expL)))/((pi_p+pi_g)*(psi*kh)) + (1-delta));


[name='Budget Constrain']

c= (A*(kh(-1)^expH)*((psi*kh(-1))^expL)) - pi_g*(psi*kh - (1-delta)*(psi*kh(-1))) - pi_p*(kh - (1-delta)*kh(-1));


[name='low-tech capital']
kl = ((pi_p*(1-varpi))/(pi_g*varpi))*kh;

[name='total capital']
k= kh + kl;

[name='investment']
i= y-c;
end;

%----------------------------------------------------------------
% 4. steadexp(y) state
%----------------------------------------------------------------
steady_state_model;

kh = (a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alpha*(1-varpi)))/(((1/beta-1+delta)*(pi_p))/(alpha*varpi)))^(1/(1-alpha));


y = a*(kh^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alpha*(1-varpi))); // output

c= (a*(kh^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alpha*(1-varpi)))) - pi_g*delta*(((pi_p*(1-varpi))/(pi_g*varpi))*kh) - pi_p*delta*kh; // consumption

kl = ((pi_p*(1-varpi))/(pi_g*varpi))*kh;


k= kh + kl;

i= y-c;

end;



%----------------------------------------------------------------
% 5. initial value
%----------------------------------------------------------------
initval;

pi_p = 1.76;
pi_g =1.47;


kh = (a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alpha*(1-varpi)))/(((1/beta-1+delta)*(pi_p))/(alpha*varpi)))^(1/(1-alpha));


y = a*(kh^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alpha*(1-varpi))); // output

c= (a*(kh^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alpha*(1-varpi)))) - pi_g*delta*(((pi_p*(1-varpi))/(pi_g*varpi))*kh) - pi_p*delta*kh; // consumption

kl = ((pi_p*(1-varpi))/(pi_g*varpi))*kh;


k= kh + kl;

i= y-c;

end;
//steady(solve_algo=4,maxit=100);
check;
//Computes the eigenvalues of the model linearized around the values specified by the last initval, endval or steady statement. Generally, the eigenvalues are only meaningful if the linearization is done around a steady state of the model. It is a device for local analysis in the neighborhoodof this steady state (Dynare Manual).


%----------------------------------------------------------------
% 6. end value (permanet shock)
%----------------------------------------------------------------
endval;

pi_p = 0.79;
pi_g =0.79;


kh = log((a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alpha*(1-varpi)))/(((1/beta-1+delta)*(pi_p))/(alpha*varpi)))^(1/(1-alpha)));


y = log(a*(exp(kh)^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh))^(alpha*(1-varpi)))); // output

c= log((a*(exp(kh)^(alpha*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh))^(alpha*(1-varpi)))) - pi_g*delta*(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh)) - pi_p*delta*exp(kh)); // consumption

kl = log(((pi_p*(1-varpi))/(pi_g*varpi))*exp(kh));


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
perfect_foresight_setup(periods=208);
perfect_foresight_solver(maxit=20,no_homotopy);


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

//dynasave(modelsimul);