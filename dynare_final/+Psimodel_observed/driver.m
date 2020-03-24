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
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'Psimodel_observed';
M_.dynare_version = '4.6.1';
oo_.dynare_version = '4.6.1';
options_.dynare_version = '4.6.1';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('Psimodel_observed.log');
M_.exo_names = cell(2,1);
M_.exo_names_tex = cell(2,1);
M_.exo_names_long = cell(2,1);
M_.exo_names(1) = {'psi_l'};
M_.exo_names_tex(1) = {'psi\_l'};
M_.exo_names_long(1) = {'psi_l'};
M_.exo_names(2) = {'psi_h'};
M_.exo_names_tex(2) = {'psi\_h'};
M_.exo_names_long(2) = {'psi_h'};
M_.endo_names = cell(5,1);
M_.endo_names_tex = cell(5,1);
M_.endo_names_long = cell(5,1);
M_.endo_names(1) = {'c'};
M_.endo_names_tex(1) = {'c'};
M_.endo_names_long(1) = {'c'};
M_.endo_names(2) = {'y'};
M_.endo_names_tex(2) = {'y'};
M_.endo_names_long(2) = {'y'};
M_.endo_names(3) = {'k'};
M_.endo_names_tex(3) = {'k'};
M_.endo_names_long(3) = {'k'};
M_.endo_names(4) = {'i'};
M_.endo_names_tex(4) = {'i'};
M_.endo_names_long(4) = {'i'};
M_.endo_names(5) = {'tfp'};
M_.endo_names_tex(5) = {'tfp'};
M_.endo_names_long(5) = {'tfp'};
M_.endo_partitions = struct();
M_.param_names = cell(9,1);
M_.param_names_tex = cell(9,1);
M_.param_names_long = cell(9,1);
M_.param_names(1) = {'varpi'};
M_.param_names_tex(1) = {'\varpi'};
M_.param_names_long(1) = {'varpi'};
M_.param_names(2) = {'betaa'};
M_.param_names_tex(2) = {'\betaa'};
M_.param_names_long(2) = {'betaa'};
M_.param_names(3) = {'al'};
M_.param_names_tex(3) = {'a_l'};
M_.param_names_long(3) = {'al'};
M_.param_names(4) = {'ah'};
M_.param_names_tex(4) = {'a_h'};
M_.param_names_long(4) = {'ah'};
M_.param_names(5) = {'sigmaa'};
M_.param_names_tex(5) = {'\sigmaa'};
M_.param_names_long(5) = {'sigmaa'};
M_.param_names(6) = {'deltaa'};
M_.param_names_tex(6) = {'\deltaa'};
M_.param_names_long(6) = {'deltaa'};
M_.param_names(7) = {'a'};
M_.param_names_tex(7) = {'(a_h^\varpi)(a_l^{1-\varpi})'};
M_.param_names_long(7) = {'a'};
M_.param_names(8) = {'alphaa'};
M_.param_names_tex(8) = {'\alphaa'};
M_.param_names_long(8) = {'alphaa'};
M_.param_names(9) = {'theta'};
M_.param_names_tex(9) = {'\dfrac{(1-\varpi)}{\varpi}'};
M_.param_names_long(9) = {'theta'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 5;
M_.param_nbr = 9;
M_.orig_endo_nbr = 5;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.orig_eq_nbr = 5;
M_.eq_nbr = 5;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 2 7;
 0 3 8;
 1 4 0;
 0 5 0;
 0 6 0;]';
M_.nstatic = 2;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 1;
M_.ndynamic   = 3;
M_.dynamic_tmp_nbr = [7; 2; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , 'Aggregate Output' ;
  2 , 'name' , 'Euler Equation' ;
  3 , 'name' , 'Budget Constrain' ;
  4 , 'name' , 'Investment' ;
  5 , 'name' , 'TFP' ;
};
M_.mapping.c.eqidx = [2 3 4 ];
M_.mapping.y.eqidx = [1 2 3 4 ];
M_.mapping.k.eqidx = [1 2 3 ];
M_.mapping.i.eqidx = [4 ];
M_.mapping.tfp.eqidx = [5 ];
M_.mapping.psi_l.eqidx = [1 5 ];
M_.mapping.psi_h.eqidx = [1 5 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [3 ];
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(5, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(9, 1);
M_.endo_trends = struct('deflator', cell(5, 1), 'log_deflator', cell(5, 1), 'growth_factor', cell(5, 1), 'log_growth_factor', cell(5, 1));
M_.NNZDerivatives = [18; -1; -1; ];
M_.static_tmp_nbr = [4; 1; 0; 0; ];
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
M_.params(7) = M_.params(4)^M_.params(1)*M_.params(3)^(1-M_.params(1));
a = M_.params(7);
M_.params(9) = (1-M_.params(1))/M_.params(1);
theta = M_.params(9);
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.exo_steady_state(2) = psi_h0;
oo_.exo_steady_state(1) = psi_l0;
oo_.steady_state(3) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8)));
oo_.steady_state(2) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*((M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8))))^M_.params(8);
oo_.steady_state(1) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*((M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8))))^M_.params(8)-(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8)))*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*M_.params(6);
oo_.steady_state(4) = oo_.steady_state(2)-oo_.steady_state(1);
oo_.steady_state(5) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1));
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
resid
oo_.dr.eigval = check(M_,options_,oo_);
%
% ENDVAL instructions
%
ys0_= oo_.steady_state;
ex0_ = oo_.exo_steady_state;
oo_.exo_steady_state(2) = psi_hF_obser;
oo_.exo_steady_state(1) = psi_lF_obser;
oo_.steady_state(3) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8)));
oo_.steady_state(2) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*((M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8))))^M_.params(8);
oo_.steady_state(1) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*((M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8))))^M_.params(8)-(M_.params(8)/(1/M_.params(2)-(1-M_.params(6))))^(1/(1-M_.params(8)))*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1))*(M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*M_.params(6);
oo_.steady_state(4) = oo_.steady_state(2)-oo_.steady_state(1);
oo_.steady_state(5) = (M_.params(4)*(1-oo_.exo_steady_state(2)))^M_.params(1)*(M_.params(3)*(1-oo_.exo_steady_state(1)))^(1-M_.params(1));
options_.periods = 180;
perfect_foresight_setup;
perfect_foresight_solver;
options_.periods = 208;
perfect_foresight_setup;
options_.no_homotopy = true;
options_.simul.maxit = 20;
perfect_foresight_solver;
resid;
save('Psimodel_observed_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('Psimodel_observed_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('Psimodel_observed_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('Psimodel_observed_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('Psimodel_observed_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('Psimodel_observed_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('Psimodel_observed_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
disp('Note: 1 warning(s) encountered in the preprocessor')
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
