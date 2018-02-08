% Darboux performances.
% author: Nicola De Renzis
% date: 29/01/2018


clear all  %#ok<CLALL>
close all

%% FIRST SECTION: ONE FIXED EIGENVALUE USING THEORETICAL EQUATION. COMPARE AGAINSt N of points



eig = [10+5*1i]; % eigenvalues
qd = [2*1i-3]; % spectral amplitude
E_NFT = 4.*sum(imag(eig));  % energy carried by discrete spectrum. In Soliton case all the energy of the signal

%Narray = 2.^(5:20);
%Narray = 2^10;

for index=1:length(Narray)
    
    % define waveform
    N = Narray(index);    % number of samples positive axis
    T1 = -20; T2 = 20; % [s]
    t = linspace(T1,T2,2*N+1); % [s]
    
    % apply INFT
    q = INFT_Darboux(eig,qd,t); % eig from theory
    if ~isrow(q)
        q=q.';
    end
    
    qd = (eig-conj(eig))*qd; % only if you want the bk
    %t1 = 1/2/imag(eig)*log(abs(qd));
    %q_th = 2*imag(eig)*sech(2*imag(eig).*(t-t1)).*exp(-1i.*(2*real(eig).*t+angle(qd)+pi));
    %t1 = 1/2/imag(eig)*log(abs(qd)/2/imag(eig));
    %q_th = 2*imag(eig)*sech(2*imag(eig).*(t-t1)).*exp(-1i.*(2*real(eig).*t+angle(qd)+pi/2));
    E_th(index) = trapz(t,abs(q_th).^2); % energy of signal;
    
    
    E(index) = trapz(t,abs(q).^2); % energy of signal;
    
    pulseError(index) = mean(abs(q_th-q).^2);
    
end
energyError = abs(E_th-E).^2;

%%{
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


%{
figure(1)
semilogy(Narray,pulseError)
title('Error comparing waveforms')
xlabel('N samples')
ylabel('MSE(q-q_{th})')

figure(2)
semilogy(Narray,energyError)
title('Error comparing energy')
xlabel('N samples')
ylabel('MSE(E-E_{th})')
%}