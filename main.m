clc;
clear all;

%% size  
data=load('ASP_Final_Data.mat');
matX=data.matX;
theta_s_noisy=data.theta_s_noisy;
theta_i_noisy=data.theta_i_noisy;

N=length(matX(:,1)); %number of sensors(N=10)
L=length(matX(1,:)); %length of time(L=2000)

t=[1:L];


%% plot noisy measurements of DOAs
figure(1)
plot(t,theta_s_noisy);
hold on
grid on
plot(t,theta_i_noisy);
hold on
legend('noisy $\theta_{s}$','noisy $\theta_{i}$','interpreter','Latex','Fontsize',20);
title('noisy measurements of DOAs','Fontsize',20);
xlabel('time');
ylabel('DOA(degree)');

%% estimate theta_s(EMD) and plot
thr=0.2;
y_theta_s=EMD(theta_s_noisy,t,thr);
theta_s_hat=y_theta_s(4,:)+y_theta_s(5,:)+y_theta_s(6,:)+y_theta_s(7,:);

% figure(2)
% plot(t,y_theta_s(1,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF1','FontSize',12);
% hold on
% figure(3)
% plot(t,y_theta_s(2,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF2','FontSize',12);
% figure(4)
% plot(t,y_theta_s(3,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF3','FontSize',12);
% figure(5)
% plot(t,y_theta_s(4,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF4','FontSize',12);
% figure(6)
% plot(t,y_theta_s(5,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF5','FontSize',12);
% figure(7)
% plot(t,y_theta_s(6,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF6','FontSize',12);
% figure(8)
% plot(t,y_theta_s(7,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF7','FontSize',12);
% figure(9)
% plot(t,theta_s_hat)
% figure(10)
% plot(t,theta_s_noisy,t,theta_s_hat)

%% estimate theta_i(EMD) and plot
y_theta_i=EMD(theta_i_noisy,t,thr);
theta_i_hat=y_theta_i(4,:)+y_theta_i(5,:);

% figure(2)
% plot(t,y_theta_i(1,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF1','FontSize',12);
% hold on
% figure(3)
% plot(t,y_theta_i(2,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF2','FontSize',12);
% figure(4)
% plot(t,y_theta_i(3,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF3','FontSize',12);
% figure(5)
% plot(t,y_theta_i(4,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF4','FontSize',12);
% figure(6)
% plot(t,y_theta_i(5,:))
% xlabel('time(sec)','FontSize',12);
% title('IMF5','FontSize',12);
% figure(7)
% plot(t,theta_i_hat)
% figure(8)
% plot(t,theta_i_noisy,t,theta_i_hat)

%% calculate steering vector
for a=1:L
    for b=1:N
        as(b,a)=exp(i*pi*(b-1)*sind(theta_s_hat(a)));
        ai(b,a)=exp(i*pi*(b-1)*sind(theta_i_hat(a)));
        ai_mi(b,a)=exp(i*(b-1)*sind(theta_i_hat(a)-55));
        ai_pl(b,a)=exp(i*(b-1)*sind(theta_i_hat(a)+55));
    end
end


%% plot estimated DOAs
figure(2)
plot(t,theta_s_hat);
hold on
grid on
plot(t,theta_i_hat);
hold on
legend('$\hat{\theta}_{s}$','$\hat{\theta}_{i}$','interpreter','Latex','Fontsize',20);
title('estimated DOAs','Fontsize',20);
xlabel('time');
ylabel('DOA(degree)');

%% array steering beamformer
y_hat_as=array_steering(matX,theta_s_hat,as);
y_hat_as_real=real(y_hat_as);
y_hat_as_imag=imag(y_hat_as);

figure(3)
subplot(2,1,1)
plot(t,y_hat_as_real);
ylim([-10,10])
title('real part of estimated source signal $\hat{s}(t)$ with array steering beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');
subplot(2,1,2)
plot(t,y_hat_as_imag);
ylim([-10,10])
title('imaginary part of estimated source signal $\hat{s}(t)$ with array steering beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');

%% uniform weighting beamformer
y_hat_uw=uniform_weighting(matX);
y_hat_uw_real=real(y_hat_uw);
y_hat_uw_imag=imag(y_hat_uw);

figure(4)
subplot(2,1,1)
plot(t,y_hat_uw_real);
ylim([-10,10])
title('real part of estimated source signal $\hat{s}(t)$ with uniform weighting beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');
subplot(2,1,2)
plot(t,y_hat_uw_imag);
ylim([-10,10])
title('imaginary part of estimated source signal $\hat{s}(t)$ with uniform weighting beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');

%% MVDR beamformer
y_hat_MVDR=MVDR(matX,theta_s_hat,as,0.1,0.99);
y_hat_MVDR_real=real(y_hat_MVDR);
y_hat_MVDR_imag=imag(y_hat_MVDR);

figure(5)
subplot(2,1,1)
plot(t,y_hat_MVDR_real);
ylim([-10,10])
title('real part of estimated source signal $\hat{s}(t)$ with MVDR beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');
subplot(2,1,2)
plot(t,y_hat_MVDR_imag);
ylim([-10,10])
title('imaginary part of estimated source signal $\hat{s}(t)$ with MVDR beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');

%% LCMV beamformer
y_hat_LCMV=LCMV(matX,theta_s_hat,theta_i_hat,as,ai,0.1,0.99);
y_hat_LCMV_real=real(y_hat_LCMV);
y_hat_LCMV_imag=imag(y_hat_LCMV);

figure(6)
subplot(2,1,1)
plot(t,y_hat_LCMV_real);
ylim([-10,10])
title('real part of estimated source signal $\hat{s}(t)$ with LCMV beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');
subplot(2,1,2)
plot(t,y_hat_LCMV_imag);
ylim([-10,10])
title('imaginary part of estimated source signal $\hat{s}(t)$ with LCMV beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');

%% my beamformer
s_t_hat=my_beamformer(matX,theta_s_hat,theta_i_hat,as,ai_mi,ai,ai_pl,0.1,0.99);
s_t_hat_real=real(s_t_hat);
s_t_hat_imag=imag(s_t_hat);

figure(7)
subplot(2,1,1)
plot(t,s_t_hat_real);
ylim([-10,10])
title('real part of estimated source signal $\hat{s}(t)$ with my beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');
subplot(2,1,2)
plot(t,s_t_hat_imag);
ylim([-10,10])
title('imaginary part of estimated source signal $\hat{s}(t)$ with my beamformer','interpreter','Latex','Fontsize',12);
xlabel('time');
ylabel('amplitude');

save("r11942135_ASPFinal_PerformanceEvaluation.mat")