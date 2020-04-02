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
           betaa $\betaa$,
           al $a_l$,

           ah $a_h$,
           sigmaa $\sigmaa$,
           deltaa $\deltaa$
           a $(a_h^\varpi)(a_l^{1-\varpi})$,
           alphaa $\alphaa$,
           theta  $\dfrac{(1-\varpi)}{\varpi}$;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
// parameters values
%Parameters RBC Cycles


 load parameterfile;
 set_param_value('betaa',betaa)
 set_param_value('varpi',varpi)
 set_param_value('al',al)
 set_param_value('ah',ah)
 set_param_value('al',al)
 set_param_value('ah',ah)
 set_param_value('sigmaa',sigmaa)
 set_param_value('deltaa',deltaa)
 set_param_value('alphaa',alphaa)

a = (ah^varpi)*(al^(1-varpi));
theta = (1-varpi)/varpi;



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;


# psi = ((pi_p*(1-varpi))/(pi_g*varpi));

# A = (ah^varpi)*(al^(1-varpi));

# expH = alphaa*varpi;

# expL = alphaa*(1-varpi);


[name='Aggregate Output']
y = A*(kh(-1)^expH)*((psi*kh(-1))^expL);


[name='Euler Equation']

c^(-sigmaa) = betaa*(c(+1)^(-sigmaa))*( (expH*(A*(kh^expH)*((psi*kh)^expL)))/(kh*(pi_p+pi_g)) + (expL*(A*(kh^expH)*((psi*kh)^expL)))/((pi_p+pi_g)*(psi*kh)) + (1-deltaa));


[name='Budget Constrain']

c= (A*(kh(-1)^expH)*((psi*kh(-1))^expL)) - pi_g*(psi*kh - (1-deltaa)*(psi*kh(-1))) - pi_p*(kh - (1-deltaa)*kh(-1));


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

kh = (a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alphaa*(1-varpi)))/(((1/betaa-1+deltaa)*(pi_p))/(alphaa*varpi)))^(1/(1-alphaa));


y = a*(kh^(alphaa*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alphaa*(1-varpi))); // output

c= (a*(kh^(alphaa*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alphaa*(1-varpi)))) - pi_g*deltaa*(((pi_p*(1-varpi))/(pi_g*varpi))*kh) - pi_p*deltaa*kh; // consumption

kl = ((pi_p*(1-varpi))/(pi_g*varpi))*kh;


k= kh + kl;

i= y-c;

end;



%----------------------------------------------------------------
% 5. initial value
%----------------------------------------------------------------
initval;

pi_p = pi_p0;
pi_g =pi_g0;


kh = (a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alphaa*(1-varpi)))/(((1/betaa-1+deltaa)*(pi_p))/(alphaa*varpi)))^(1/(1-alphaa));


y = a*(kh^(alphaa*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alphaa*(1-varpi))); // output

c= (a*(kh^(alphaa*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alphaa*(1-varpi)))) - pi_g*deltaa*(((pi_p*(1-varpi))/(pi_g*varpi))*kh) - pi_p*deltaa*kh; // consumption

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

pi_p = pi_pF_simul;
pi_g =pi_gF_simul;

kh = (a*((((pi_p*(1-varpi))/(pi_g*varpi)))^(alphaa*(1-varpi)))/(((1/betaa-1+deltaa)*(pi_p))/(alphaa*varpi)))^(1/(1-alphaa));


y = a*(kh^(alphaa*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alphaa*(1-varpi))); // output

c= (a*(kh^(alphaa*varpi))*((((pi_p*(1-varpi))/(pi_g*varpi))*kh)^(alphaa*(1-varpi)))) - pi_g*deltaa*(((pi_p*(1-varpi))/(pi_g*varpi))*kh) - pi_p*deltaa*kh; // consumption

kl = ((pi_p*(1-varpi))/(pi_g*varpi))*kh;


k= kh + kl;

i= y-c;
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