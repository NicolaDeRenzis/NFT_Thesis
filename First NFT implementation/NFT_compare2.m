% main script for comparing the NFT routine with the analytical results
% this script compare a hyperbolic secant pulse (Satsuma-Yajima signals)
clear all %#ok<CLALL>
close all

% define waveform
N = 2^12;    % number of samples
A = 1;    % amplitude
T1 = -7; T2 = 7; % [s]
t = linspace(T1,T2,2*N+1); % [s]
%t = (-N:N).*1e-2;
q = A.*sech(t);
dt = t(2)-t(1); % time step [s]
T = t(N)-t(1);

% analytical results

y = sym('y');

lambda_th = linspace(-40,40,1e5);

a = @(x) (gamma(-1i.*x + 0.5)).^2 ./ (gamma(-1i.*x + 0.5 + A).*gamma(-1i.*x + 0.5 - A));
b = @(x) 1i.*sin(pi.*A)./cosh(pi.*x);
di=1e-3;
imagAxis = [0:di:10].*1i;

qc_theory = b(lambda_th)./a(lambda_th);

a_disc = 1;
for i=1:2:100
    if A-i/2>0
        r = (i-1)/2+1;
        lambdaj_th(r) = 1i*(A-i/2);
        bk(r) = 1i*(-1)^(r-1); % VALID ONLY IF A IS INTEGER
        a_disc = a_disc.*((imagAxis - lambdaj_th(r))./(imagAxis - conj(lambdaj_th(r))));
    end
end
da = gradient(a_disc);
[~,ind] = find(imagAxis==lambdaj_th.');
qd_th = bk./da(ind);



% NFT numerical - Forward
[qc,lambda,qd,lambdaj] = NFT_search(q,t,'fw');

%% Plots
figure(1)
plot(lambda_th,abs(qc_theory));
hold on
plot(lambda, abs(qc),'.')
hold off
xlabel('\lambda')
ylabel({'$|\hat{q}(\lambda)|$'},'Interpreter','latex')
grid on
legend('Analytical','FW NFT')
title('Continuous NFT Spectrum')

figure(2)
plot(lambdaj_th,'o','markersize',10)
hold on
plot(lambdaj,'o','markersize',10)
xlabel('Re')
ylabel('Im')
grid on
legend('Analytical','FW NFT')
title('Discrete NFT Spectrum')
hold off

figure(3)
plot(t,q)
title('pulse')
xlabel('t [s]')
ylabel('q(t)')