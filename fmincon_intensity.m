%Function to optmise the geometry for a zone plate helium microscope.

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
% figure;
% x_min=min(xData);
% x_max=max(xData);
% x_plot=linspace(x_min,x_max,10000);
% y_plot=10.^(fitresult(x_plot));
% x_plot=10.^(x_plot);
% 
% h = plot( x_plot,y_plot,'r');
% hold on
% plot(10.^(xData), 10.^(yData),'bx')
% %legend( h, 'S_log vs. P_d_log', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel('P_0d_{noz} / torr cm')
% ylabel('Terminal parallel speed ratio')
% set(gca,'Xscale','log','Yscale','log')
% xlim([1 400])
% ylim([5 500])
% %grid on

%% Create model for intensity and beam width
reset(symengine)
syms S S_log P_d_log P d_noz B n lambda_wav f delta_r lambda sigma tau
m=4.002602*1.660539040e-27;
k_B=1.38064852e-23;
T=300;

P_d_log=log10(P*d_noz);

%Create model for speed ratio using the result from the fit
S_log=fitresult.a*P_d_log+fitresult.b+fitresult.c./(1+exp(-fitresult.d*(P_d_log-fitresult.mu)));
S=10.^S_log;


%Brightness can also be looked at
Br=P/S;
Br_diff=diff(Br,P);
Br_h=matlabFunction(subs(Br));
Br_diff_h=matlabFunction(subs(Br_diff));
S_h=matlabFunction(subs(S));
conv_mbar_factor=1.33322;

%Plotting brightness
% d_plot=10e-4; %Nozzle diameter in cm
% P_plot1=0:10:130e3*(10e-4/d_plot);
% P_plot2=130e3*(10e-4/d_plot):10:2e6;
% 
% 
% figure;h1=plot(P_plot1*conv_mbar_factor*1e-3,0.18*conv_mbar_factor*100*Br_h(P_plot1,10e-4)/sqrt(m*k_B*T),'b');
% hold on
% plot(P_plot2*conv_mbar_factor*1e-3,0.18*conv_mbar_factor*100*Br_h(P_plot2,10e-4)/sqrt(m*k_B*T),'b--')
% xlabel('Nozzle pressure/bar')
% ylabel('Source brightness/s^{-1}sr^{-1}m^{-2}')
% set(gca,'XScale','log')
% xlim([1 1500])
% 
% d_plot=5e-4;
% P_plot1=0:10:130e3*(10e-4/d_plot);
% P_plot2=130e3*(10e-4/d_plot):10:2e6;
% h2=plot(P_plot1*conv_mbar_factor*1e-3,0.18*conv_mbar_factor*100*Br_h(P_plot1,d_plot)/sqrt(m*k_B*T),'r');
% plot(P_plot2*conv_mbar_factor*1e-3,0.18*conv_mbar_factor*100*Br_h(P_plot2,d_plot)/sqrt(m*k_B*T),'r--')
% 
% d_plot=20e-4;
% P_plot1=0:10:130e3*(10e-4/d_plot);
% P_plot2=130e3*(10e-4/d_plot):10:2e6;
% h3=plot(P_plot1*conv_mbar_factor*1e-3,0.18*conv_mbar_factor*100*Br_h(P_plot1,d_plot)/sqrt(m*k_B*T),'m');
% plot(P_plot2*conv_mbar_factor*1e-3,0.18*conv_mbar_factor*100*Br_h(P_plot2,d_plot)/sqrt(m*k_B*T),'m--')
% 
% legend([h2 h1 h3],{'5\mum nozzle','10\mum nozzle','20\mum nozzle'},'Location','northwest')


%Zone plate radius
rN=lambda_wav*f/(2*delta_r);

mult=1e20;
%Negative intensity in beam. So that when we minimise we find the max 
I=-(n*P*rN^2*B^2/S) *mult;%* (0.18*pi^2/sqrt(m*k_B*T))


%Beam width condition
g=-(sigma-sqrt((0.42*delta_r)^2+(B*f/sqrt(3))^2+(rN/(S*sqrt(2)*tau))^2));



%% Find constrained maximum using fmincon
%Set all the constants
d_noz=10e-4;        %Nozzle diameter in cm
delta_r=25e-9;      %Minimum feature size of zone plate
f=1e-3;             %Working Distance
n=0.05;             %Efficiency of the zone plate
tau=3.582;          %Chromatic aberration corretion term
lambda_wav=5.63e-11;%Helium wavelength
sigma_loop=[delta_r*1.23:0.5e-9:80e-9,80e-9:5e-9:2000e-9];
N_loop=length(sigma_loop);
x=zeros(N_loop,2);
I_out=zeros(N_loop,1);

%Loop over every resolution
for k=1:N_loop
    
    %Set the target resolution for this loop
    sigma=sigma_loop(k);
    
    %Write matlab functions of the intensity and constraints
    matlabFunction(subs(I),'File','I_func','Vars',[B P]);
    matlabFunction(subs(g),'File','g_func','Vars',[B P]);
    
    %Set the functions to be used by optimiser. Inputs need manipulating to
    %be accepted. Also scale the parameters to make it easier for
    %optimiser.
    I_fun=@I_func_mb;
    nonlcon=@g_func_mb;
    
    %Set options
    options=optimoptions('fmincon','Display','iter');
    %options.ConstraintTolerance=1e-;
    
    %Find the constrained minimum
    [x(k,:),I_out(k)]=fmincon(I_fun,[sigma/(1e-5*f),0.4],[0,1],1.75,[],[],[0,0],[],nonlcon,options);
end

%Change back to absolute intensity
I_out=I_out* (0.18*pi^2/sqrt(m*k_B*T))*conv_mbar_factor*100/mult; 

x(:,2)=x(:,2)*1e5*conv_mbar_factor/1e3; %Change back to bar

S_h=matlabFunction(subs(S));
Br_diff_h=matlabFunction(subs(Br_diff));

%Find the pressure which gives the optimum brightness
P_o=fzero(Br_diff_h,4e4);

%Calculate analytic approximation
gamma=(P_o*conv_mbar_factor*100)/(S_h(P_o))* (0.18*pi^2/(4*sqrt(m*k_B*T)));
F_th=3*(gamma*n*lambda_wav^2/delta_r^2)*(sigma_loop.^2-(0.42*delta_r)^2);

%Plot difference between analyic and numerical solutions
%figure;plot(sigma_loop,(F_th'-(-I_out))./F_th')


