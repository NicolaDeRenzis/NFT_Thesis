% Darboux with noisy eigs.
% author: Nicola De Renzis
% date: 05/02/2018


clear all  %#ok<CLALL>
close all

%% FIRST SECTION: ONE FIXED EIGENVALUE USING THEORETICAL EQUATION. USE NOISY EIGS

% define eigenvalue
eig_clean = [1+0.5*1i -1+1i*0.5 1.5*1i]; % eigenvalues
qd_clean = [1-1i 1+1i 2*1i]; % spectral amplitude

% add noise
var_i = 1e-2;    % variance in imaginary axis
var_r = 1e-2;    % variance in real axis
n = 1e3;         % number of noisy samples
noise = randn(n,length(eig_clean))*sqrt(var_r) + 1i.*randn(n,length(eig_clean)).*sqrt(var_i);
noisy_eig = repmat(eig_clean,n,1) + noise;
noisy_eig(imag(noisy_eig)<=0) = conj(noisy_eig(imag(noisy_eig)<=0));
noisy_qd = repmat(qd_clean,n,1)+noise;

%%{ debug figure
figure(1)
plot(real(noisy_eig),imag(noisy_eig),'b.')
grid on
%}

% define waveform parameters
N = 2^12; % number of samplez per signal
T1 = -12; T2 = 12; % [s]
t = linspace(T1,T2,2*N+1); % [s]

q = zeros(2*N+1,length(n));
q_th = zeros(2*N+1,length(n));
q_clean = INFT_Darboux(eig_clean,qd_clean,t);
for index=1:n
    disp(index)
    eig = noisy_eig(index,:);
    qd = noisy_qd(index,:);
    
    % apply INFT
    q(:,index) = INFT_Darboux(eig,qd,t);
    
    %qd = (eig-conj(eig))*qd; % only if you want the bk
    %t1 = 1/2/imag(eig)*log(abs(qd));
    %q_th = 2*imag(eig)*sech(2*imag(eig).*(t-t1)).*exp(-1i.*(2*real(eig).*t+angle(qd)+pi));
    %t1 = 1/2/imag(eig)*log(abs(qd)/2/imag(eig));
    %q(:,index) = 2*imag(eig)*sech(2*imag(eig).*(t-t1)).*exp(-1i.*(2*real(eig).*t+angle(qd)+pi/2));
    
    % useless stuff
    %E_th(index) = trapz(t,abs(q_th).^2); % energy of signal
    %E(index) = trapz(t,abs(q).^2); % energy of signal
    %E_NFT = 4.*sum(imag(eig));  % energy carried by discrete spectrum. In Soliton case all the energy of the signal
    
end


% figures
close all
figure(2)
subplot(2,2,1:2)
plot(t,abs(q).^2)
ylabel('|q(t)|^2')
xlabel('t')
subplot(2,2,3:4)
plot(t,angle(q)./pi)
ylabel('phase/\pi')
xlabel('t')


figure(3)
subplot(2,2,[1,3])
plot(0,0,'b',0,0,'k',t,real(q),'b',t,real(q_clean),'k') % first two only for legend purposes
ylabel('q(t)')
xlabel('t')
ylim([-4,2])
xlim([T1,T2])
legend('Real part','Clean real part','Location','Southeast')
grid on
subplot(2,2,[2,4])
plot(0,0,'r',0,0,'k',t,imag(q),'r',t,imag(q_clean),'k') % first two only for legend purposes
ylabel('q(t)')
xlabel('t')
ylim([-4,2])
xlim([T1,T2])
legend('Imag part','Clean imag part','Location','Southeast')
grid on

figure(4)
plot(t,real(mean(q,2)),'b',t,imag(mean(q,2)),'r--',t,real(q_clean),'k',t,imag(q_clean),'k--')
ylabel('q(t)')
xlabel('t')
legend('Real part','Imag part','Clean signal')
grid on

%{
figure(4)
plot(t,real(q_th),'b',t,imag(q_th),'b--')
hold on
plot(t,real(q),'r',t,imag(q),'r--')
%plot(t,abs(q_th).^2,'g--',t,abs(q).^2,'k-.')
hold off
ylabel('|q(t)|^2')
xlabel('t [s]')
legend('Analytic - Real','Analytic - Imag','Darboux - Real','Darboux - Imag')

figure(5)
subplot(2,1,2)
plot(t,real(q_th)-real(q),'b');title('Real part')
ylabel('Error')
subplot(2,1,1)
plot(t,imag(q_th)-imag(q),'r');title('Imag part')
ylabel('Error')
xlabel('time')
%}

