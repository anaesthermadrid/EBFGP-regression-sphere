clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONAL GAUSSIAN PROCESS (FGP) REGRESSION BASED ON TIME-ADAPTIVE EMPIRICAL BAYES (EB) IN THE SPHERE
%
% OBJECTIVE: 
% This code provides the implementation of FGP under time-adaptive EB approach 
% from FGP prior, whose time-varying covariance kernel is constructed from
% Matérn  spatiotemporal covariance family (in the Gneiting class)
% restricted to the sphere in R3, displaying
% Short-Range Dependence (SRD) and  possible fractality. 
% 
%--------------------------------------------------------------------------
% VARIANCE AND BIAS ANALYSIS OF THE POSTERIOR FUNCTIONAL PREDICTOR:
%
% OUTPUTS (Numerical results and visualization) 
%---------------------------------------------------------------------------------------------------------
% Posterior functional parameters of the infinite-dimensional predictive Gaussian probability distribution
%---------------------------------------------------------------------------------------------------------
% - Functional series defining the posterior mean
% - Posterior covariance operators at the observed times
%--------------------------------------------------------------------------
% Posterior functional variance and bias analysis
%--------------------------------------------------------------------------
% -Empirical mean quadratic functional errors
% -Posterior time linear correlation
% -Functional bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 1: SPATIAL MESH AND ANGLE NODES CONFIGURATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- STEP 1.1: RESOLUTION AND TRUNCATION PARAMETERS STEP 1.1: RESOLUTION AND TRUNCATION PARAMETERS (UNDER SPHERICAL SPARSE OBSERVATION OF FUNCTIONAL DATA) --- ---

T = 100;           % Functional sample size (number of observed times)
TR = fix(log(T));  % Number of eigenspaces (truncation parameter value given by logarithmic rule).
                   % Optionally one can consider the power function truncation scheme, TR=T^(beta), 0<beta <1
N = 50;            % Number of spherical locations 
                 
M = 100;           % Number of nodes at each component 

%% --- STEP 1.2: SAMPLING ANGULAR FREQUENCY ---

% 1.2.1. Polar and azimuthal nodes
rvals = linspace(-1, 1, N)';             % Grid in the Legendre support  [-1,1]
elevation = asin(rvals);                 % Polar angle nodes in (-pi/2, pi/2)
anguloP = elevation + pi/2 * ones(N, 1); % Polar angle nodes in [0, pi] 
azimuth = linspace(0, 2*pi, N)';         % Azimuthal angle nodes in [0, 2pi]

% 1.2.2. Two-dimensional meshgrid
[Theta, Phi] = meshgrid(anguloP, azimuth); % Theta and Phi (N x N)

% 1.2.3. Transform angle meshgrid outputs into Cartesian coordinates
X = sin(Theta) .* cos(Phi);
Y = sin(Theta) .* sin(Phi);
Z = cos(Theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 2: PRIOR COVARIANCE HYPERPARAMETER SAMPLING AND MATÉRN SPATIOTEMPORAL COVARIANCE  SUBFAMILY RESTRICTED TO SPHERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- STEP 2.1: INFORMATIVE HYPERPARAMETER PRIOR DISTRIBUTION PROPOSAL AND SAMPLING ---

% 2.1.1 Degenerate parameter scenarios
a = 1;
c = 1; 
nabla = 1;

% 2.1.2  Samples drawn from hyperparameter priors

% FGP memory covariance hyperparameters
A1 = 8; B1 = 2; beta = betarnd(A1,B1,M,1);
A2 = 11; B2 = 5; alpha = betarnd(A2,B2,M,1);

% FGP local spatial regularity covariance hyperparameters
mu = 5;
sigma = 0.1; 
nu = normrnd(mu,sigma,M,1);

% Additive unstructured observation noise variance hyperparameter (contamination parameter)    
mu = 2;
sigma = 0.01; 
sigma = normrnd(mu,sigma,M,1);
 
%% --- STEP 2.2: RESTRICTION OF MATÉRN COVARIANCE FUNCTION TO THE SPHERE (TRUE HYPERPARAMETER VECTOR VALUE) ---
 
% 2.2.1 Spatial distance matrix (chordal distance on the sphere)   
C1 = zeros(N, N);
for i = 1:N
    C1(i,i) = 1;
    for j = 1:i-1
        C1(i,j) = 2 * sin((anguloP(i) - anguloP(j)) / 2);
        C1(j,i) = C1(i,j);
    end
end

% 2.2.2 Evaluation of positive function psi with complete monotone derivatives
TAU = linspace(0, T, T)'; % Time lag vector
theta0 = unidrnd(M); % Index for true covariance hyperparameter vector
psi=(ones(T,1)+a*TAU.^(2*alpha (theta0,1) )).^(beta (theta0,1)); % Function psi
 
% 2.2.3. Spatiotemporal covariance kernel restricted to sphere (Subfamily A3 Gneiting Class)
psinv=(psi).^(-1);

% Array dimension
CEB=zeros(N, N, T);

for i=1:N            
    for j=1:N
        scnorm=C1(i,j)*psinv;
        CEB(i,j,:)=psinv.*(nabla^(2)*(2^(nu(theta0,1)-1)*gamma (nu(theta0,1)))^(-1)*(c*scnorm).^(nu (theta0,1))).*besselk(nu(theta0,1),c*scnorm);       
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 3: TIME-VARYING ANGULAR SPECTRUM AND KL EXPANSION VIA TRUNCATED SPHERICAL HARMONIC TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% --- STEP 3.1: ELEMENTS OF THE  TRUNCATED SPHERICAL HARMONIC BASIS ---
% Generation of the real-valued zonal spherical harmonics using Legendre polynomials.

% Array dimension
RR1 = zeros(N, TR);

for n = 1:TR
    % Normalizing factor (zonal case m=0)
    % Considering the second index m=0 in all Laplace-Beltrami eigenspaces
    % since  working with an isotropic model in the sphere
    factor_norm = sqrt((2*n+1) / (4*pi)); 
    x_leg = cos(anguloP);  % Argument for Legendre polynomial evaluation 
    
    % Legendre polynomial evaluation and corresponding values of spherical harmonics
    Pn = legendreP(n, x_leg); %  Nx1 dimension
    Yl0 = factor_norm * Pn;   % Real-valued spherical harmonic (zonal)
    RR1(:,n) = Yl0;           % Store in basis matrix (N x TR)
end

%% --- STEP 3.2: GRAPHICAL REPRESENTATION OF ELEMENTS OF TRUNCATED SPHERICAL HARMONIC BASIS FUNCTIONS ---

zzb1 = 1;
figure('Name', 'Truncated spherical harmonic basis');
sgtitle('Truncated spherical harmonic basis')
for n = 1:TR     
    SHF = RR1(:,n);
    SHF_mesh = repmat(SHF.', N, 1);           
    subplot(ceil(TR/2), 2, zzb1)
    surf(X, Y, Z, SHF_mesh);
    colormap(jet(256)); colorbar('vert'); shading interp
    title(['$n=$ ', num2str(n)], 'Interpreter', 'latex')
    hold on 
    zzb1 = zzb1 + 1;
end
     
%% --- STEP 3.3: TIME-VARYING ANGULAR SPECTRUM ---

% Array dimension
TEIGC = zeros(T, TR);

% Time-varying FGP angular spectrum 
for n=1:TR
    for t=1:T   
        CC= CEB(:,:,t);         
        TEIGC(t,n)=RR1(:,n)'*CC*RR1(:,n);        
    end
end
   
% Array dimension 
TASDEB = zeros(TR, TR, T);

for t = 1:T  
    TASDEB(:,:,t) = diag(TEIGC(t,:)');
end

%% --- STEP 3.4: TIME-VARYING KARHUNEN-LOÈVE-BASED SIMULATION OF FGP PRIOR---
 
R = 200; % Number of replicates of FGP prior involved in the implementation of the time-adaptive EB approach 
         % (minimum threshold for implementation of Monte Carlo numerical integration in time-adaptive EB)
r0 = unidrnd(R); % The choice made in the generated replicates for graphical representation
                 % of FGP, additive observation noise, observation process, functional bias

% Array dimension
KLRCOREB = zeros(TR, T, R); % FGP  time-varying random coefficients (TR eigenspaces x T times x R replicates)
TORIGEB = zeros(TR, N, T);  % (TR eigenspaces x N spherical locations x T times) input argument dimension
YORIGEB = zeros(N, T, R);   % Inverse truncated spherical harmonic transform (N spherical locations x T times x R replicates)

for r = 1:R
    for t = 1:T
        KLRCOREB(:,t,r)= (abs(TEIGC(t,:))').^(1/2).*randn(TR,1); 
        for j = 1:N
            TORIGEB(:,j,t) = KLRCOREB(:,t).*RR1(j,:)';
        end
    end   
    STORIGEB = sum(TORIGEB);
    YORIGEB(:,:,r) = STORIGEB(1,:,:);
end
 
%% --- STEP 3.5: VISUALIZATION OF FGP PRIOR  ---
 
zzb2 = 1;
DELTA = T/8 + 1; % Sampling times in graphical representation
nrows = ceil(length(1:fix(DELTA):T)/2);
figure('Name', 'FGP prior');
sgtitle('FGP prior')
for t = 1:fix(DELTA):T   
    Yl02 = YORIGEB(:,t,r0);
    Yl02_mesh = repmat(Yl02.', N, 1);   
    subplot(nrows,2,zzb2)
    surf(X, Y, Z,Yl02_mesh);
    colormap(jet(256)); colorbar('vert'); shading interp
    title(['$t=$ ', num2str(t)], 'Interpreter', 'latex')
    hold on 
    zzb2 = zzb2 + 1;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 4: ADDITIVE OBSERVATION NOISE, AND OBSERVATION PROCESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- STEP 4.1: ADDITIVE GAUSSIAN OBSERVATION NOISE GENERATION ---
 
% Generation of time-varying random coefficients based on equality in probability distribution
GWNOBS = sigma(theta0) * randn(TR, T, R); % theta0 provides the true value which can be selected from the support of the corresponding hyperparameter prior

% Array dimension
GWN2 = zeros(TR, N, T);  
SEGWN = zeros(N, T, R);

% Inverse truncated spherical harmonic transform
for r = 1:R
    for t = 1:T
        for i = 1:N        
            GWN2(:,i,t) = GWNOBS(:,t,r) .* RR1(i,:)';
        end
        SEWN = sum(GWN2, 1);  
        SEGWN(:,t,r) = SEWN(1,:,t);
    end
end
  
%% --- STEP 4.2: VISUALIZATION OF OBSERVATION NOISE ---

zzb3 = 1;
figure('Name', 'Truncated generation of additive observation noise');
sgtitle('Truncated generation of additive observation noise')
for t = 1:fix(DELTA):T        
    AON = SEGWN(:,t,r0); 
    GAON = repmat(AON.', N, 1);             
    subplot(nrows,2,zzb3)
    surf(X,Y,Z, GAON);
    colormap(jet(256)); colorbar('vert'); shading interp
    title(['$t= $ ', num2str(t)], 'Interpreter', 'latex')
    hold on 
    zzb3 = zzb3 + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 5: TIME ADAPTIVE EMPIRICAL BAYES (ML-II)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- STEP 5.1: EVALUATION OF SPATIOTEMPORAL COVARIANCE FUNCTIONS AT THE PRIOR HYPERPARAMETER SAMPLES GENERATED ---
% Evaluation of psi function
 
for theta = 1:M
    psiEB(:,theta) = (ones(T,1) + a * TAU.^(2 * alpha(theta,1))).^(beta(theta,1));
end

% Array dimension
CEBB = zeros(N, N, T, M);

% Evaluation of  Matérn subfamily example of spatiotemporal covariance
% function, restricted to the sphere
for theta = 1:M 
    psinv_k = (psiEB(:,theta)).^(-1);
    for i = 1:N            
        for j = 1:N
            scnorm = C1(i,j) * psinv_k;
             CEBB(i,j,:,theta) = psinv_k.*(nabla^(2)*(2^(nu(theta,1)-1)*gamma (nu(theta,1)))^(-1)*(c*scnorm).^(nu (theta,1))).*besselk(nu(theta,1),c*scnorm);    
        end
    end
end
  
%--------------------------------------------------------------------------
% Evaluation of time-varying angular spectrum of FGP 
%--------------------------------------------------------------------------

% Array dimension
TEIGC_EB = zeros(T, TR, M);

% Computation of time-varying angular spectrum at the generated hyperparameter prior samples
for theta = 1:M
    for n = 1:TR
        for t = 1:T
            CC = CEBB(:,:,t,theta);
            TEIGC_EB(t,n,theta) = RR1(:,n)' * CC * RR1(:,n); 
        end
    end
end

%% --- STEP 5.2: TIME-VARYING ANGULAR SPECTRUM OF OBSERVATION PROCESS ---
 
%  Array dimension
NCM = zeros(TR, TR, T, M);  
TASD = zeros(TR, TR, T, M);  
 
% Time-varying diagonal matrixes associated with the truncated spectral decomposition
% of the involved covariance operator families of the FGP and the observation noise
for theta = 1:M
    for t = 1:T
        ANCM = diag(sigma(theta)^(2)*ones(TR,1)); 
        ATASD = diag(TEIGC_EB(t,:,theta)');           
        NCM(:,:,t,theta) = ANCM; 
        TASD(:,:,t,theta) = ATASD;
    end
end

% Computation of the time-varying diagonal matrix associated with the
% diagonal purely spectral decomposition of the covariance operator family
% of the functional Gaussian observation process
PPSOBS = NCM + TASD; 

 
%% --- STEP 5.3: TRUNCATED FREDHOLM DETERMINANT OF COVARIANCE OPERATORS OF FGP AND OBSERVATION PROCESS ---
 

TR2 = 120;       % Truncation of Fredholm determinants in the number of powers of time-varying 
                 % covariance operators involved in their computations for FGP and observation process
PERT = 0.01;     % Regularization parameter for numerical stability defining suitable neighborhood
EPSILON = PERT * ones(T,M);

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Time-varying truncated Fredholm determinant of the covariance operator
% family of the FGP
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


% Array dimension
TRAZTV = zeros(T, M);

% Time-varying truncated trace norm of the covariance operator family of
% the FGP
for theta = 1:M
    for t = 1:T
        TRAZTV(t,theta) = sum(diag(TASD(:,:,t,theta))); 
    end
end

for theta = 1:M
    ITRAZA2(:,theta) = TRAZTV(:,theta).^(-1)';
    OMEGA2(:,theta) = ITRAZA2(:,theta) - EPSILON(:,theta);
end

%----------------------------------------------------------------------------------------------------------------------------------
% Time-varying truncated Fredholm determinant of the covariance operator family of the FGP
%----------------------------------------------------------------------------------------------------------------------------------

for theta = 1:M
    for k = 1:TR2
        v1B = TRAZTV(:,theta); 
        v2B = OMEGA2(:,theta);
        FD02(k,:) = exp(-(v1B.^(k).*v2B.^(k))/k)';
    end
    TVFDOR(:,theta) = (sum(FD02))'; % Truncated Time-Varying Fredholm Determinant of FGP
end

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
% Time-varying truncated Fredholm determinant of the covariance operator family of the functional observation process
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

% Array dimension
TRAZA = zeros(T, M);

for theta = 1:M
    for t = 1:T 
        TRAZA(t,theta) = sum(diag(PPSOBS(:,:,t,theta))); 
    end
end 

 
for theta = 1:M
    ITRAZA(:,theta) = TRAZA(:,theta).^(-1)';
    OMEGA(:,theta) = ITRAZA(:,theta) - EPSILON(:,theta);
end

 
for theta = 1:M
    for k = 1:TR2
        v1 = TRAZA(:,theta); 
        v2 = OMEGA(:,theta);
        FD0(k,:) = exp(-(v1.^(k).*v2.^(k))/k)';
    end
    TVFD(:,theta) = (sum(FD0))'; % Truncated Time-Varying Fredholm Determinant of observation process
end
  
%% ---STEP 5.4: COMPUTATION OF ML-II PARAMETER ESTIMATE---

% Computation of truncated log-likelihood of observation process depending on FGP and observation noise intensity hyperparameter
for r = 1:R
    for theta = 1:M
        for t = 1:T          
            IRKHSP(t,theta,r) = (KLRCOREB(:,t,r) + GWNOBS(:,t,r))' * diag((diag(PPSOBS(:,:,t,theta))).^(-1)) * (KLRCOREB(:,t,r) + GWNOBS(:,t,r));          
            LL(t,theta,r) = (-1/2) * (log(TVFD(t,theta))) + IRKHSP(t,theta,r);
        end
    end
end

% Computation of truncated log-likelihood of FGP depending on its covariance function hyperparameter vector 
for theta = 1:M
    for r = 1:R
        for t = 1:T  
            IRKHSTHETA(t,theta,r) = KLRCOREB(:,t,r)' * diag((diag(TASD(:,:,t,theta))).^(-1)) * KLRCOREB(:,t,r);             
            LLTHETA(t,theta,r) = (-1/2) * (log(TVFDOR(t,theta)) + IRKHSTHETA(t,theta,r));
        end
    end
end

% Computation, via Monte Carlo numerical integration, of truncated likelihood of observation process 
% depending on the FGP hyperparameter covariance, and observation noise intensity  
MLD = mean(exp(LL).*exp(LLTHETA), 3);

% Computation of ML-II parameter vector estimate 

% Array dimension
MLII = zeros(T,1);

% Maximum of  truncated likelihood of observation process 
% depending on the FGP hyperparameter covariance, and observation noise intensity
for t=1:T
    [~, idx] = max(MLD(t,:)');
    MLII(t) = idx;  
end 
MLII = MLII';
 
%% --- STEP 5.5: CONDITIONAL TIME-VARYING SPECTRUM GIVEN ML-II ESTIMATOR VALUES ---

% Array dimension
TASDEBE = zeros(TR,TR,T);   
PPSOBSEB = zeros(TR,TR,T);

% Evaluation of FGP  and observation process time-varying spectrumm
% conditionally to the time-adaptive ML-II computed estimates
for t = 1:T
    theta = MLII(t); 
    TASDEBE(:,:,t) = TASD(:,:,t, theta);  
    PPSOBSEB(:,:,t) = PPSOBS(:,:,t, theta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STAGE 6: TRUNCATED INFINITE-DIMENSIONAL PREDICTIVE PROBABILITY DISTRIBUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- STEP 6.1: TRUNCATED POSTERIOR FUNCTIONAL MEAN (SEE THEOREM 1 OF RUIZ-MEDINA et al., 2026) --- 

for r = 1:R   
    for t = 1:T
        PROYFPOSTMEANEB(:,t,r)=TASDEBE(:,:,t)*(inv(PPSOBSEB(:,:,t)))*(KLRCOREB(:,t,r)+GWNOBS(:,t,r));   
        FPOSTMEAN2EB(r,t,:,:)=RR1*diag(PROYFPOSTMEANEB(:,t,r));   
    end
end
 
% Reconstruction via inverse spherical harmonic transform (N=500 is the minimum 
% threshold to obtaining reasonable spherical discretization errors)

FFEB = sum(FPOSTMEAN2EB,4); 

%--------------------------
% Graphical representation 
%--------------------------

% Array dimension
FPOSTMEAN2VFEB = zeros(N,T,R); 

for r = 1:R
    for t = 1:T
        for i = 1:N           
            FPOSTMEAN2VFEB(i,t,r) = FFEB(r,t,i);
        end
    end
end

zzb4 = 1;
sgtitle('Posterior mean predictor')
figure('Name', 'Posterior Mean Predictor');
for t = 1:fix(DELTA):T       
    Yl03 = FPOSTMEAN2VFEB(:,t,r0);
    Yl03_mesh = repmat(Yl03.', N, 1);       
    subplot(nrows,2,zzb4)
    surf(X, Y, Z,Yl03_mesh);
    colormap(jet(256)); colorbar('vert'); shading interp
    title(['$t=$ ', num2str(t)], 'Interpreter', 'latex')
    hold on 
    zzb4 = zzb4 + 1;
end
 
%% --- STEP 6.2: TIME-VARYING PROJECTED EMQE ---

% Truncated projected EMQE over time (l2 identification, see Theorem 1 in Ruiz-Medina et al., 2026)
EMQE= mean((PROYFPOSTMEANEB - KLRCOREB).^(2), 3);   

%------------------------------------------------------
% Graphical representation 
%------------------------------------------------------
zzb5 = 1;
figure('Name', 'Projected EMQE over time');
sgtitle('Projected EMQE over time')
for n = 1:TR
    EMQEEIG = EMQE(n,:)';
    subplot(ceil(TR/2),2,zzb5)
    plot(1:T, EMQEEIG, 'go',...
        'MarkerSize', 3, 'MarkerFaceColor', [0.5,0.5,0.5]);
    title(['$n=$ ', num2str(n)], 'Interpreter', 'latex')
    xlabel('Time'); ylabel('EMQE');
    hold on; 
    zzb5 = zzb5 + 1;
end

%Summary statistics providing time-average of above EMQE
ATEMQE = mean(EMQE, 2); %  

%------------------------------------------------------
% Graphical representation
%------------------------------------------------------

figure('Name', 'Averaged in time of EMQEs');
sgtitle('Averaged in time of EMQEs')
plot(1:TR, ATEMQE, '--go',...
    'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', [0.5,0.5,0.5]);
xlabel('Eigenspace'); ylabel('ATEMQE');

 
%% --- STEP 6.3: POSTERIOR APPROXIMATION OF TIME LINEAR CORRELATION OF FGP MODEL ---

% Array dimension
TV2 = zeros(TR,T);

for n = 1:TR
    for t = 1:T        
        TV2(n,t) = var(KLRCOREB(n,t,:));  
    end
end

%------------------------------------------------------
% Graphical representation 
%------------------------------------------------------
 
zzb6 = 1;
figure('Name', 'Posterior time linear correlation of FGP');
sgtitle('Posterior time linear correlation of FGP')
for n = 1:TR
    TVEIG = (TV2(n,6:T))'; % Theoretical  linear correlation in time for functional GP (t=6 being minimum threshold to remove temporal edge effect)
    TVEIG2 = (TV2(n,6:T) - EMQE(n,6:T))'; % Posterior linear correlation in time for FGP  (t=6 being minimum threshold to remove temporal edge effect)
    subplot(ceil(TR/2),2,zzb6)
    plot(6:T, TVEIG, 'bo', 6:T, TVEIG2, 'ro');
    legend('Theoretical', 'Posterior');
    title(['$n=$ ', num2str(n)], 'Interpreter', 'latex');
    xlabel('Time'); 
    hold on; 
    zzb6 = zzb6 + 1;
end   
      
%% --- STEP 6.4: TIME-VARYING SPHERICAL FUNCTIONAL BIAS AVERAGE BASED ON  R-REPLICATES OF THE POSTERIOR PREDICTOR  ---

% Array dimension
BIAS = zeros(TR, T, R);
BBB  = zeros(TR, N, T);
SPHBIAS = zeros(N,T);

%--------------------------------------------------------------------------------------------------------------
%Computation of functional spherical bias series averaged over the R replicates
%--------------------------------------------------------------------------------------------------------------

for r = 1:R   
    for t = 1:T
        BIAS(:,t,r) =  KLRCOREB(:,t,r)-TASDEBE(:,:,t)*(pinv(PPSOBSEB(:,:,t)))*(KLRCOREB(:,t,r));         
    end
end 

% Average over R replicates
BB = (mean(BIAS,3)); 

for t = 1:T
    for j = 1:N         
        BBB(:,j,t) = BB(:,t).*RR1(j,:)';
    end
end

BBBB = sum(BBB,1);  
SPHBIAS(:,:) = BBBB(1,:,:);

%------------------------------------------
% Graphical representation 
%------------------------------------------

zzb7 = 1;
figure('Name', 'Time-varying spherical functional averaged bias');
sgtitle('Time-varying spherical functional averaged bias')
DELTA = fix(T/8)+1;
for t = 1:DELTA:T   
    Yl04 = SPHBIAS(:,t); 
    Yl04_mesh = repmat(Yl04.', N, 1);   
    subplot(nrows,2,zzb7)
    surf(X,Y,Z,Yl04_mesh);
    colormap(jet(256)); colorbar('vert'); shading interp
    title(['$t=$ ',num2str(t)], 'Interpreter','latex')
    hold on 
    zzb7 = zzb7+1;
 end
 