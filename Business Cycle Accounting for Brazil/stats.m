% New stats
% This program computes correlations and standard deviations of the model

% Compute Output in log-deviations

ythat = log(yt) - log(mean(yt));

% Compute model predictions in log-deviations

yt_at = (log(yt_at) - log(mean(yt_at)))';
yt_tht = (log(yt_tht) - log(mean(yt_tht)))';
yt_tkt = (log(yt_tkt) - log(mean(yt_tkt)))';
yt_tbt = (log(yt_tbt) - log(mean(yt_tbt)))';


%Run HP filter on series to separate trend from cycle


v_obs = [ythat yt_at yt_tht yt_tkt yt_tbt at tht tkt tbt];
v_filter = zeros(size(v_obs));

for i = 1:9
v = v_obs(:,i);
hp_filter;
v_filter(:,i) = vhp;
end

% 1) First Describe the properties of the Wedges

disp('This is the vector of relative Standard Deviations of Wedges (1995-2018');
sd_wedges= [std(v_filter)];
sd_wedges = sd_wedges/sd_wedges(1);
SD =sd_wedges(6:9);
T = array2table(SD, 'VariableNames',{'TFP', 'Labor_Wedge', 'Capital_Wedge', 'Bond_Wedge'});
T