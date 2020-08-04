
%%TVP-VAR Heterodedástico (baseado no codigo do Koop)
%Variveis: Taxa de Varição do PIB e FBCF (% PIB) sobre o PIB e variação do PIB
%Codigo para o artigo da disciplina de Econometria II PPGEco/UFSC 2017
%Denise Manfredini


randn('state',sum(100*clock));
rand('twister',sum(100*clock));
%% ----------------------------------DADOS----------------------------------------

load ydata.dat
load yearlab.dat

% you have to x13 the series
Y = ydata;
%yearlab =datenum(2003,1:208,1)';

% Dimensões de X e Y
t=size(Y,1); % t is the time-series observations of Y
M=size(Y,2); % M is the dimensionality of Y

% Defasagens
p = 1; % defasagem do VAR
numa = M*(M-1)/2; % Number of lower triangular elements of A_t (other than 0's and 1's)
%%  ===================================| VAR EQUATION |==============================
% Gera a matriz de regressores X_t = [1 y_t-1 y_t-2 ... y_t-k] para t=1:T
ylag = mlag2(Y,p); % Y is [T x M]. ylag is [T x (Mp)]
ylag = ylag(p+1:t,:);

K = M + p*(M^2); % K é o número de elementos no vetor de estados
% Cria a matriz Z_t.
Z = zeros((t-p)*M,K);
for i = 1:t-p
    ztemp = eye(M);
    for j = 1:p        
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];  %#ok<AGROW>
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end

% Redefine y
y = Y(p+1:t,:)';
yearlab = yearlab(p+1:t);

% Tamanho da série temporal
t=size(y,2);  

%% ----------------------------PRELIMINARES---------------------------------
% Preliminares para o Gibbs
nrep = round(0.80*500000);  % Amostragens
nburn = round(0.20*nrep);   % burn-in
it_print = round(0.05*nrep);

% Escolha dos anos da FIR
FIR1 = 1963;
FIR2 = 1979;
FIR3 = 2015;

%========= PRIORS NAO INFORMATIVA:

 A_OLS = zeros(numa,1);
 B_OLS = zeros(K,1);
 VA_OLS = eye(numa);
 VB_OLS = eye(K);
 sigma_OLS = [-9; 0; 0]; %ones tinha dados certo

sizeW = M; % Size of matrix W
sizeS = 1:M; % Size of matrix S

%-------- PRIORES

% B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prmean = B_OLS;
B_0_prvar = 4*VB_OLS;
% A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prmean = A_OLS;
A_0_prvar = 4*VA_OLS;
% log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prmean = sigma_OLS;
sigma_prvar = eye(M);

% Note that for IW distribution I keep the _prmean/_prvar notation....
% Q is the covariance of B(t), S is the covariance of A(t) and W is the
% covariance of (log) SIGMA(t)
% Q ~ IW(k2_Q*size(subsample)*Var(B_OLS),size(subsample))
Q_prmean = ((k_Q)^2)*VB_OLS;
Q_prvar = 10;
% W ~ IG(k2_W*(1+dimension(W))*I_n,(1+dimension(W)))
W_prmean = ((k_W)^2)*ones(M,1);
W_prvar = 2;
% S ~ IW(k2_S*(1+dimension(S)*Var(A_OLS),(1+dimension(S)))
S_prmean = cell(M-1,1);
S_prvar = zeros(M-1,1);
ind = 1;
for ii = 2:M
   
    S_prmean{ii-1} = ((k_S)^2)*(1 + sizeS(ii-1))*VA_OLS(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind);
    S_prvar(ii-1) = 1 + sizeS(ii-1);
    ind = ind + ii;
end

%% ========= MARIZES DE INICIALIZACAO:

consQ = 0.0001;
consS = 0.0001;
consH = 0.01;
consW = 0.0001;
Ht = kron(ones(t,1),consH*eye(M));   
Htchol = kron(ones(t,1),sqrt(consH)*eye(M)); 
Qdraw = consQ*eye(K);  
Sdraw = consS*eye(numa); 
Sblockdraw = cell(M-1,1); 
ijc = 1;
for jj=2:M
    Sblockdraw{jj-1} = Sdraw(((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc);
    ijc = ijc + jj;
end
Wdraw = consW*ones(M,1);  
Btdraw = zeros(K,t);     
Atdraw = zeros(numa,t);  
Sigtdraw = zeros(t,M);  
sigt = kron(ones(t,1),0.01*eye(M));   
statedraw = 5*ones(t,M);       
                              
Zs = kron(ones(t,1),eye(M));

% Para guardar os resultados
Bt_postmean = zeros(K,t);   
At_postmean = zeros(numa,t); 
Sigt_postmean = zeros(t,M);  
Qmean = zeros(K,K);          % Q
Smean = zeros(numa,numa);    % S
Wmean = zeros(M,1);          % W
Bt_save=NaN(K,t,nrep);

sigmean = zeros(t,M);   
cormean = zeros(t,numa); 
sig2mo = zeros(t,M);     
cor2mo = zeros(t,numa); 

%========= IMPULSO RESPOSTA:

istore = 1;
if istore == 1
    nhor = 15;  % Impulse response horizon
    imp75 = zeros(nrep,M,nhor);
    imp81 = zeros(nrep,M,nhor);
    imp96 = zeros(nrep,M,nhor);
    bigj = zeros(M,M*p);
    bigj(1:M,1:M) = eye(M);
end


%% ====================================== INICIO DO AMOSTRADOR SAMPLING ========================================

tic; % This is just a timer
disp('Number of iterations');


for irep = 1:nrep + nburn    % Inicio do Gibbs
  
    if mod(irep,it_print) == 0
        disp(irep);toc;
    end
    % -----------------------------------------------------------------------------------------
    %   STEP I: Amostra B de p(B|y,A,Sigma,V) (Coeficientes dos estados)
    % -----------------------------------------------------------------------------------------


    draw_beta
    
     Bt_save(:,:,irep)=Btdraw;
    
    %-------------------------------------------------------------------------------------------
    %   STEP II: Amostra A(t) de p(At|y,B,Sigma,V) (Coeficientes dos Estados)
    %-------------------------------------------------------------------------------------------
    
    
    draw_alpha
    
    
    %------------------------------------------------------------------------------------------
    %   STEP III: Amostra a matriz de covariância do VAR log-SIGMA(t)
    %------------------------------------------------------------------------------------------
    
    draw_sigma
    
    
    % Create the VAR covariance matrix H(t). It holds that:
    %           A(t) x H(t) x A(t)' = SIGMA(t) x SIGMA(t) '
    Ht = zeros(M*t,M);
    Htsd = zeros(M*t,M);
    for i = 1:t
        inva = inv(capAt((i-1)*M+1:i*M,:));
        stem = diag(sigt(i,:));
        Hsd = inva*stem;
        Hdraw = Hsd*Hsd';
        Ht((i-1)*M+1:i*M,:) = Hdraw;  % H(t)
        Htsd((i-1)*M+1:i*M,:) = Hsd;  % Cholesky DE H(t)
    end
    
    %----------------------------IMPULSO RESPOSTAS E RESULTADOS PARA DESCARTE-----------------
    if irep > nburn      
        
        Bt_save(:,:,nrep) = Btdraw;
        Bt_postmean = Bt_postmean + Btdraw;   % regression coefficients B(t)
        At_postmean = At_postmean + Atdraw;   % lower triangular matrix A(t)
        Sigt_postmean = Sigt_postmean + Sigtdraw;  % diagonal std matrix SIGMA(t)
        Qmean = Qmean + Qdraw;     % covariance matrix Q of B(t)
        ikc = 1;
        for kk = 2:M
            Sdraw(((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc)=Sblockdraw{kk-1};
            ikc = ikc + kk;
        end
        Smean = Smean + Sdraw;    % covariance matrix S of A(t)
        Wmean = Wmean + Wdraw;    % covariance matrix W of SIGMA(t)
        % Correlacoes e variancias variantes no tempo
        stemp6 = zeros(M,1);
        stemp5 = [];
        stemp7 = [];
        for i = 1:t
            stemp8 = corrvc(Ht((i-1)*M+1:i*M,:));
            stemp7a = [];
            ic = 1;
            for j = 1:M
                if j>1
                    stemp7a = [stemp7a ; stemp8(j,1:ic)']; 
                    ic = ic+1;
                end
                stemp6(j,1) = sqrt(Ht((i-1)*M+j,j));
            end
            stemp5 = [stemp5 ; stemp6']; 
            stemp7 = [stemp7 ; stemp7a']; 
        end
        sigmean = sigmean + stemp5; % diagonal da matriz de covariância do VAR
        cormean =cormean + stemp7;  % elementos fora da diagonal dos elementos da matriz de covariância do VAR
        sig2mo = sig2mo + stemp5.^2;
        cor2mo = cor2mo + stemp7.^2;
         
        if istore==1
            
            
            IRA_tvp
            
            
        end 
    end 
end 
clc;
toc; 
%============================= TERMINA O AMOSTRADOR DE GIBBS ==================================
Bt_postmean = Bt_postmean./nrep;  % Posterior mean of B(t) (VAR regression coeff.)
At_postmean = At_postmean./nrep;  % Posterior mean of A(t) (VAR covariances)
Sigt_postmean = Sigt_postmean./nrep;  % Posterior mean of SIGMA(t) (VAR variances)
Qmean = Qmean./nrep;   % Posterior mean of Q (covariance of B(t))
Smean = Smean./nrep;   % Posterior mean of S (covariance of A(t))
Wmean = Wmean./nrep;   % Posterior mean of W (covariance of SIGMA(t))



sigmean = sigmean./nrep;
cormean = cormean./nrep;
sig2mo = sig2mo./nrep;
cor2mo = cor2mo./nrep;


