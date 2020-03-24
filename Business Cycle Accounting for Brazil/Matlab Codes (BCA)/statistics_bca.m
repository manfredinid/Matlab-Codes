% This program computes correlations and standard deviations of the model

% Compute Output in log-deviations

ythat = log(yt) - log(mean(yt));

% Compute model predictions in log-deviations

yt_at = (log(yt_at) - log(mean(yt_at)))';
yt_tht = (log(yt_tht) - log(mean(yt_tht)))';
yt_tkt = (log(yt_tkt) - log(mean(yt_tkt)))';
yt_tbt = (log(yt_tbt) - log(mean(yt_tbt)))';

% Filter all series 

v=ythat;
hp_filter;
ythat=vhp;

v=yt_at;
hp_filter;
yt_at=vhp;

v=yt_tht;
hp_filter;
yt_tht=vhp;

v=yt_tkt;
hp_filter;
yt_tkt=vhp;

v=yt_tbt;
hp_filter;
yt_tbt=vhp;

v=at;
hp_filter;
at=vhp;

v=tht;
hp_filter;
tht=vhp;

v=tkt;
hp_filter;
tkt=vhp;

v=tbt;
hp_filter;
tbt=vhp;

% 1) First Describe the properties of the Wedges

disp('This is the vector of relative Standard Deviations of Wedges');
disp('TFP, Labor Wedge, Capital Wedge, Bond Wedge (1991-2006)');

sd_wedges=[std(at)/std(ythat);std(tht)/std(ythat);std(tkt)/std(ythat);std(tbt)/std(ythat)]


disp('This is the vector of Correlation +- 2 lags');
disp('TFP, Labor Wedge, Capital Wedge, Bond Wedge (1991-2006)');

a=corrcoef(at((1+2):T),ythat(1:(T-2)));
b=corrcoef(tht((1+2):T),ythat(1:(T-2)));
c=corrcoef(tkt((1+2):T),ythat(1:(T-2)));
d=corrcoef(tbt((1+2):T),ythat(1:(T-2)));

corr_neg2=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(at((1+1):T),ythat(1:(T-1)));
b=corrcoef(tht((1+1):T),ythat(1:(T-1)));
c=corrcoef(tkt((1+1):T),ythat(1:(T-1)));
d=corrcoef(tbt((1+1):T),ythat(1:(T-1)));

corr_neg1=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(at((1):T),ythat(1:(T)));
b=corrcoef(tht((1):T),ythat(1:(T)));
c=corrcoef(tkt((1):T),ythat(1:(T)));
d=corrcoef(tbt((1):T),ythat(1:(T)));

corr_0=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(at(1:(T-1)),ythat((1+1):T));
b=corrcoef(tht(1:(T-1)),ythat((1+1):T));
c=corrcoef(tkt(1:(T-1)),ythat((1+1):T));
d=corrcoef(tbt(1:(T-1)),ythat((1+1):T));

corr_1=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(at(1:(T-2)),ythat((1+2):T));
b=corrcoef(tht(1:(T-2)),ythat((1+2):T));
c=corrcoef(tkt(1:(T-2)),ythat((1+2):T));
d=corrcoef(tbt(1:(T-2)),ythat((1+2):T));

corr_2=[a(1,2); b(1,2); c(1,2); d(1,2)];

corr_wedges=[corr_neg2, corr_neg1, corr_0, corr_1, corr_2]

% 2) Then Describe the properties of the models prediction

disp('This is the vector of relative Standard Deviations of Model Prediction');
disp('TFP, Labor Wedge, Capital Wedge, Bond Wedge (1991-2006)');

sd_model=[std(yt_at)/std(ythat);std(yt_tht)/std(ythat);std(yt_tkt)/std(ythat);std(yt_tbt)/std(ythat)]


disp('This is the vector of Correlation of Model Prediction +- 2 lags');
disp('TFP, Labor Wedge, Capital Wedge, Bond Wedge (1991-2006)');

a=corrcoef(yt_at((1+2):T),ythat(1:(T-2)));
b=corrcoef(yt_tht((1+2):T),ythat(1:(T-2)));
c=corrcoef(yt_tkt((1+2):T),ythat(1:(T-2)));
d=corrcoef(yt_tbt((1+2):T),ythat(1:(T-2)));

corr_neg2=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(yt_at((1+1):T),ythat(1:(T-1)));
b=corrcoef(yt_tht((1+1):T),ythat(1:(T-1)));
c=corrcoef(yt_tkt((1+1):T),ythat(1:(T-1)));
d=corrcoef(yt_tbt((1+1):T),ythat(1:(T-1)));

corr_neg1=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(yt_at((1):T),ythat(1:(T)));
b=corrcoef(yt_tht((1):T),ythat(1:(T)));
c=corrcoef(yt_tkt((1):T),ythat(1:(T)));
d=corrcoef(yt_tbt((1):T),ythat(1:(T)));

corr_0=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(yt_at(1:(T-1)),ythat((1+1):T));
b=corrcoef(yt_tht(1:(T-1)),ythat((1+1):T));
c=corrcoef(yt_tkt(1:(T-1)),ythat((1+1):T));
d=corrcoef(yt_tbt(1:(T-1)),ythat((1+1):T));

corr_1=[a(1,2); b(1,2); c(1,2); d(1,2)];

a=corrcoef(yt_at(1:(T-2)),ythat((1+2):T));
b=corrcoef(yt_tht(1:(T-2)),ythat((1+2):T));
c=corrcoef(yt_tkt(1:(T-2)),ythat((1+2):T));
d=corrcoef(yt_tbt(1:(T-2)),ythat((1+2):T));

corr_2=[a(1,2); b(1,2); c(1,2); d(1,2)];

corr_model=[corr_neg2, corr_neg1, corr_0, corr_1, corr_2]

disp('Coefficients / Standard Errors / T-Statistics');

x1

x1se

ts=x1./x1se;
ts