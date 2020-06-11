%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'model_simul';
M_.dynare_version = '4.5.7';
oo_.dynare_version = '4.5.7';
options_.dynare_version = '4.5.7';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('model_simul.log');
M_.exo_names = 'pi_g';
M_.exo_names_tex = 'pi\_g';
M_.exo_names_long = 'pi_g';
M_.exo_names = char(M_.exo_names, 'pi_p');
M_.exo_names_tex = char(M_.exo_names_tex, 'pi\_p');
M_.exo_names_long = char(M_.exo_names_long, 'pi_p');
M_.endo_names = 'kh';
M_.endo_names_tex = 'kh';
M_.endo_names_long = 'kh';
M_.endo_names = char(M_.endo_names, 'kl');
M_.endo_names_tex = char(M_.endo_names_tex, 'kl');
M_.endo_names_long = char(M_.endo_names_long, 'kl');
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'y');
M_.endo_names_tex = char(M_.endo_names_tex, 'y');
M_.endo_names_long = char(M_.endo_names_long, 'y');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.endo_partitions = struct();
M_.param_names = 'varpi';
M_.param_names_tex = '\varpi';
M_.param_names_long = 'varpi';
M_.param_names = char(M_.param_names, 'betaa');
M_.param_names_tex = char(M_.param_names_tex, '\betaa');
M_.param_names_long = char(M_.param_names_long, 'betaa');
M_.param_names = char(M_.param_names, 'al');
M_.param_names_tex = char(M_.param_names_tex, 'a_l');
M_.param_names_long = char(M_.param_names_long, 'al');
M_.param_names = char(M_.param_names, 'ah');
M_.param_names_tex = char(M_.param_names_tex, 'a_h');
M_.param_names_long = char(M_.param_names_long, 'ah');
M_.param_names = char(M_.param_names, 'sigmaa');
M_.param_names_tex = char(M_.param_names_tex, '\sigmaa');
M_.param_names_long = char(M_.param_names_long, 'sigmaa');
M_.param_names = char(M_.param_names, 'deltaa');
M_.param_names_tex = char(M_.param_names_tex, '\deltaa');
M_.param_names_long = char(M_.param_names_long, 'deltaa');
M_.param_names = char(M_.param_names, 'a');
M_.param_names_tex = char(M_.param_names_tex, '(a_h^\varpi)(a_l^{1-\varpi})');
M_.param_names_long = char(M_.param_names_long, 'a');
M_.param_names = char(M_.param_names, 'alphaa');
M_.param_names_tex = char(M_.param_names_tex, '\alphaa');
M_.param_names_long = char(M_.param_names_long, 'alphaa');
M_.param_names = char(M_.param_names, 'theta');
M_.param_names_tex = char(M_.param_names_tex, '\dfrac{(1-\varpi)}{\varpi}');
M_.param_names_long = char(M_.param_names_long, 'theta');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 6;
M_.param_nbr = 9;
M_.orig_endo_nbr = 6;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('model_simul_static');
erase_compiled_function('model_simul_dynamic');
M_.orig_eq_nbr = 6;
M_.eq_nbr = 6;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 1 2 0;
 0 3 0;
 0 4 8;
 0 5 0;
 0 6 0;
 0 7 0;]';
M_.nstatic = 4;
M_.nfwrd   = 1;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 1;
M_.nspred   = 1;
M_.ndynamic   = 2;
M_.equations_tags = {
  1 , 'name' , 'Aggregate Output' ;
  2 , 'name' , 'Euler Equation' ;
  3 , 'name' , 'Budget Constrain' ;
  4 , 'name' , 'low-tech capital' ;
  5 , 'name' , 'total capital' ;
  6 , 'name' , 'investment' ;
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(6, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(9, 1);
M_.NNZDerivatives = [24; 0; -1];
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
M_.params( 7 ) = M_.params(4)^M_.params(1)*M_.params(3)^(1-M_.params(1));
a = M_.params( 7 );
M_.params( 9 ) = (1-M_.params(1))/M_.params(1);
theta = M_.params( 9 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.exo_steady_state( 2 ) = pi_p0;
oo_.exo_steady_state( 1 ) = pi_g0;
oo_.steady_state( 1 ) = (M_.params(7)*((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1)))^((1-M_.params(1))*M_.params(8))/(oo_.exo_steady_state(2)*(1/M_.params(2)-1+M_.params(6))/(M_.params(1)*M_.params(8))))^(1/(1-M_.params(8)));
oo_.steady_state( 4 ) = M_.params(7)*oo_.steady_state(1)^(M_.params(1)*M_.params(8))*((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*oo_.steady_state(1))^((1-M_.params(1))*M_.params(8));
oo_.steady_state( 3 ) = M_.params(7)*oo_.steady_state(1)^(M_.params(1)*M_.params(8))*((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*oo_.steady_state(1))^((1-M_.params(1))*M_.params(8))-(1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*oo_.steady_state(1)*oo_.exo_steady_state(1)*M_.params(6)-oo_.steady_state(1)*oo_.exo_steady_state(2)*M_.params(6);
oo_.steady_state( 2 ) = (1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*oo_.steady_state(1);
oo_.steady_state( 5 ) = oo_.steady_state(1)+oo_.steady_state(2);
oo_.steady_state( 6 ) = oo_.steady_state(4)-oo_.steady_state(3);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
oo_.dr.eigval = check(M_,options_,oo_);
%
% ENDVAL instructions
%
ys0_= oo_.steady_state;
ex0_ = oo_.exo_steady_state;
oo_.exo_steady_state( 2 ) = pi_pF_simul;
oo_.exo_steady_state( 1 ) = pi_gF_simul;
oo_.steady_state( 1 ) = log((M_.params(7)*((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1)))^((1-M_.params(1))*M_.params(8))/(oo_.exo_steady_state(2)*(1/M_.params(2)-1+M_.params(6))/(M_.params(1)*M_.params(8))))^(1/(1-M_.params(8))));
oo_.steady_state( 4 ) = log(M_.params(7)*exp(oo_.steady_state(1))^(M_.params(1)*M_.params(8))*((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*exp(oo_.steady_state(1)))^((1-M_.params(1))*M_.params(8)));
oo_.steady_state( 3 ) = log(M_.params(7)*exp(oo_.steady_state(1))^(M_.params(1)*M_.params(8))*((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*exp(oo_.steady_state(1)))^((1-M_.params(1))*M_.params(8))-oo_.exo_steady_state(1)*M_.params(6)*(1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*exp(oo_.steady_state(1))-oo_.exo_steady_state(2)*M_.params(6)*exp(oo_.steady_state(1)));
oo_.steady_state( 2 ) = log((1-M_.params(1))*oo_.exo_steady_state(2)/(M_.params(1)*oo_.exo_steady_state(1))*exp(oo_.steady_state(1)));
oo_.steady_state( 6 ) = log(exp(oo_.steady_state(4))-exp(oo_.steady_state(3)));
options_.periods = 180;
perfect_foresight_setup;
perfect_foresight_solver;
options_.periods = 208;
perfect_foresight_setup;
options_.no_homotopy = 1;
options_.simul.maxit = 20;
perfect_foresight_solver;
resid;
save('model_simul_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('model_simul_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('model_simul_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('model_simul_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('model_simul_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('model_simul_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('model_simul_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
