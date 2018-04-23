%Function to optimise the geometry for both the pinhole and zone plate
%geometries and compare them.

%Perform numeric optimisation of zone plate
fmincon_intensity
%Perform numeric optimisation of the pinhole
fmincon_intensity_pinhole


%Plot the intensities from both cases
figure;plot(sigma_loop*1e9,-I_out,'k')
hold on
plot(sigma_loop_pinhole*1e9,-I_out_pinhole,'r')
%plot(sigma_loop,F_th,'m')
xlabel('\sigma /nm')
ylabel('Optimised Beam Flux /s^{-1}')
legend('Zone Plate','Pinhole','location','southeast')
set(gca,'YScale','log')
set(gca,'XScale','log')
ylim([1e6 1e12])
xlim([35 2000])
xticks([50 100 200 300 500 1000])


%Plot the SNR for both cases with and without zero order stop.
trans_prob=0.005;
det_prob=0.001;
open_frac=0.374;
tot_prob=trans_prob*det_prob;
figure;plot(sigma_loop*1e9,sqrt(-I_out*tot_prob),'k')
hold on
plot(sigma_loop*1e9,sqrt(n/open_frac)*sqrt(-I_out*tot_prob),'b')
plot(sigma_loop_pinhole*1e9,sqrt(-I_out_pinhole*tot_prob),'r')
xlabel('\sigma /nm')
ylabel('Simulated Signal to Noise Ratio')
legend('Zone Plate with zero order stop','Zone Plate without zero order stop','Pinhole','location','northwest')
set(gca,'YScale','log')
set(gca,'XScale','log')
ylim([1 1e3])
xlim([35 2000])
xticks([50 100 200 300 500 1000])

%Angular source size plot
figure;plot(sigma_loop*1e9,x(:,1)*1e-5)
hold on
plot(sigma_loop_pinhole*1e9,x_pinhole(:,1)*1e-5)
%plot(sigma_loop_pinhole*1e9,(sqrt(3)/f)*sqrt((sigma_loop_pinhole.^2/2)-((0.42*lambda_wav*f)^2./(2*sigma_loop_pinhole.^2))))
xlabel('Beam Width \sigma/nm')
ylabel('Optimum Angular Source Size \beta /m^{-1}')
legend('Zone Plate','Pinhole','location','NorthWest')
axis tight


%Pinhole diameter plot
figure;plot(sigma_loop_pinhole*1e9,x_pinhole(:,3)*1e-6*1e9)
xlabel('Beam Width \sigma/nm')
ylabel('Optimum Pinhole Diameter\sl d\rm/nm')
axis tight

%Plot the optimal pressure
% figure;plot(sigma_loop*1e9,x(:,2))
% hold on
% plot(sigma_loop*1e9,x(:,2),'m')
% plot([1e-5, 2000],[(P_o*conv_mbar_factor*100)/1e5, (P_o*conv_mbar_factor*100)/1e5],'r--')
% xlabel('\sigma/nm')
% ylabel('Optimum Pressure P_{0}/bar')
% xlim([30 1000])
% ylim([55 100])
% xticks([50 100 200 300 500 1000])
% set(gca,'XScale','log')
% legend('Zone plate with standard ab.','Zone plate with corrected ab.','Pinhole')