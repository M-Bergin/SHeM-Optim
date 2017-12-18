%Function to optmise the geometry for a pinhole helium microscope.

%% Fit a model to speed ratio data

%Load in the data from grabit
load('Speed_ratio_data.mat')
S_log_raw=S_data(:,2);
P_d_log_raw=S_data(:,1);
S_raw=10.^S_log_raw;
P_d_raw=10.^P_d_log_raw;


%Fit the data

[xData, yData] = prepareCurveData( P_d_log_raw, S_log_raw );

% Set up fittype and options.
ft = fittype( 'a*x+b+c/(1+exp(-d*(x-mu)))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [-Inf -Inf 0 0 0];
% opts.StartPoint = [0.5 0.7 0.8 5 1.96];

opts.Lower = [-Inf -Inf -Inf 0 -Inf];
opts.StartPoint = [0.5 0.0675 0.756944648237544 0.706046088019609 0.27692298496089];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'S_log vs. P_d_log', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel P_d_log
% ylabel S_log
% grid on

%% Create model for intensity and beam width
reset(symengine)
syms S S_log P_d_log P d_noz B lambda_wav f lambda sigma d

m=4.002602*1.660539040e-27;
k_B=1.38064852e-23;
T=300;
conv_mbar_factor=1.33322;

P_d_log=log10(P*d_noz);

%Create model for speed ratio using the result from the fit
S_log=fitresult.a*P_d_log+fitresult.b+fitresult.c./(1+exp(-fitresult.d*(P_d_log-fitresult.mu)));
S=10.^S_log;


%Brightness can also be looked at
Br=P/S;
Br_diff=diff(Br,P);
Br_h=matlabFunction(subs(Br));
Br_diff_h=matlabFunction(subs(Br_diff));



%Negative intensity in beam. So that when we minimise we find the max 
I=-B^2*d^2*P/(S)*1e20; %* (0.18*pi^2/sqrt(m*k_B*T)*conv_mbar_factor*100)
%Beam width condition
g=sigma-sqrt((d/(2*sqrt(3)))^2+(B*f/sqrt(3))^2+(0.42*lambda_wav*f/d)^2);



%% Find constrained maximum using fmincon
%Set all the constants
d_noz=10e-4; %Nozzle diameter in cm
f=1e-3;
lambda_wav=5e-11;

%Set the target resolution here. Note for a pinhole the best possible
%resolution is given by sqrt(1.22*f*lambda_wav)
sigma_loop_pinhole=[sqrt(0.4205*f*lambda_wav)/3^(0.25):0.5e-9:150e-9,150e-9:5e-9:2e-6];
N_loop=length(sigma_loop_pinhole);
x_pinhole=zeros(N_loop,3);
I_out_pinhole=zeros(N_loop,1);

%Loop over every resolution
for k=1:N_loop
    
    %Set the target resolution for this loop
    sigma=sigma_loop_pinhole(k);
    
    %Write matlab functions of the intensity and constraints
    matlabFunction(subs(I),'File','I_func_pinhole','Vars',[B P d]);
    matlabFunction(subs(g),'File','g_func_pinhole','Vars',[B P d]);
    
    %Set the functions to be used by optimiser. Inputs need manipulating to
    %be accepted. Also scale the parameters to make it easier for
    %optimiser.
    I_fun=@I_func_pinhole_mb;
    nonlcon=@g_func_pinhole_mb;
    
    %Set options
    options=optimoptions('fmincon','Display','iter');
    options.MaxFunctionEvaluations=5000;
    options.ConstraintTolerance=1e-8;
    options.OptimalityTolerance=1e-8;
    %options.StepTolerance=1e-12;
    
    %Find the constrained minimum
    [x_pinhole(k,:),I_out_pinhole(k)]=fmincon(I_fun,[sigma/(1e-5*f),0.4,1.2*sigma*1e6],[0,1,0],1.25,[],[],[0,0,0],[],nonlcon,options);
end

x_pinhole(:,2)=x_pinhole(:,2)*1e5*conv_mbar_factor/1e3; %Convert back to bar

%Convert back to absolute intensity
I_out_pinhole=I_out_pinhole*1e-20* (0.18*pi^2/(4*sqrt(m*k_B*T))*conv_mbar_factor*100);

S_h=matlabFunction(subs(S));
Br_diff_h=matlabFunction(subs(Br_diff));

%Find the pressure at the optimum brightness
P_o=fzero(Br_diff_h,4e4);

%Calculate the analytic optimum flux
gamma=(P_o*conv_mbar_factor*100)/(S_h(P_o))* (0.18*pi^2/(4*sqrt(m*k_B*T)));
F_th_pinhole=3*(gamma)*(3*sigma_loop_pinhole.^4/f^2-(0.42*lambda_wav)^2);

%Plot the difference between numeric and analytic answers
figure;plot(sigma_loop_pinhole,(F_th_pinhole'-(-I_out_pinhole))./F_th_pinhole')
