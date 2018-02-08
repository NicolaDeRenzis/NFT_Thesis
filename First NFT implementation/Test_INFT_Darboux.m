% Tester for INFT_Darboux
% author: Nicola De Renzis
% date: 23/01/2018

clear all %#ok<CLALL>
close all

%% NFT of signal (eig pos from theory, eig sepctral ampl from simulation)

% define waveform
N = 2^15;    % number of samples positive axis
A = 2;    % amplitude
T1 = -16; T2 = 16; % [s]
t = linspace(T1,T2,2*N+1); % [s]
q_th = A.*sech(t).*exp(+1i*0);
h = t(2)-t(1); % time step [s]
E_th = trapz(t,abs(q_th).^2); % energy of signal; better approx of integral.

% analytical results
lambda_th = linspace(-40,40,1e5);

a = @(x) (gamma(-1i.*x + 0.5)).^2 ./ (gamma(-1i.*x + 0.5 + A).*gamma(-1i.*x + 0.5 - A));
b = @(x) 1i.*sin(pi.*A)./cosh(pi.*x);
imagAxis = [0:1e-3:10].*1i;

% continuous spectrum
qc_theory = b(lambda_th)./a(lambda_th);
E_cont = 1/pi .* trapz(lambda_th,log(1+abs(qc_theory).^2)); % energy carried by continuous spectrum

% discrete spectrum
a_disc=1;
for i=1:2:100
    if A-i/2>0
        r = (i-1)/2+1;
        lambdaj(r) = 1i*(A-i/2);
        qd_th(r) = 1i*(-1)^(r-1);
        a_disc = a_disc.*((imagAxis - lambdaj(r))./(imagAxis - conj(lambdaj(r))));
    end
end
da = gradient(a_disc);
[~,ind] = find(imagAxis==lambdaj.');
%qd_th = qd_th./da(ind);
E_disc = 4.*sum(imag(lambdaj));  % energy carried by discrete spectrum

E_NFT = E_cont + E_disc; % total energy from NFT domain
[qc,~,qd_th,lambdaj] = NFT_search(q_th,t,'fw');

%{
% nice visualization of surface a
figure(1000)
[X,Y] = meshgrid(-2:0.01:2,[0:0.01:2.1].*1i);
provaX = a([X]);
provaY = a([Y]);
mesh(X,imag(Y),real(provaY).*real(provaX))
hold on
mesh(X,imag(Y),real(provaY).*real(provaX).*0,'EdgeColor','k')
hold off
%}
%% Plots
figure(1)
plot(lambda_th,abs(qc_theory));
xlabel('\lambda')
ylabel({'$|\hat{q}(\lambda)|$'},'Interpreter','latex')
grid on
legend('Analytical')
title('Continuous NFT Spectrum')

figure(2)
plot(lambdaj_th,'o','markersize',10)
hold on
plot(lambdaj,'o','markersize',10)
hold off
xlabel('Re')
ylabel('Im')
grid on
legend('Analytical','Algorithm')
title('Discrete NFT Spectrum')
hold off

figure(3)
plot(t,q_th)
title('original pulse')
xlabel('t [s]')
ylabel('q(t)')


%% Darboux
% set eignevalues and norming constant


%q = INFT_Darboux(lambdaj_th,qd_th,t); % eig from theory
q_eig_sim = INFT_Darboux(lambdaj,qd_th,t); % eig from NFT

E = trapz(t,abs(q_eig_sim).^2); % energy of signal;

% figures
figure(4)

plot(t,real(q_th),'b',t,imag(q_th),'b--')
hold on
%plot(t,real(q),'r',t,imag(q),'r--')
plot(t,real(q_eig_sim),'g',t,imag(q_eig_sim),'g--')
%plot(t,abs(q).^2,'r--')
plot(t,abs(q_eig_sim).^2,'g-x')
plot(t,abs(q_th).^2,'b-x')
hold off
ylabel('|q(t)|^2')
xlabel('t [s]')
legend('Analytic - Real','Analytic - Imag','Darboux from algorithm - Real','Darboux from algorithm - Imag')

