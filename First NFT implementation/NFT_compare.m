% main script for comparing the NFT routine with the analytical results
% this script compare rectangular pulse 
clear all %#ok<CLALL>
close all

% define waveform
N = 2^9;    % number of samples
A = 12;      % amplitude
q = A.*ones(1,N);
T1 = 0; T2 = 1; % [s]
t = linspace(T1,T2,N); % [s]
dt = t(2)-t(1); % time step [s]
T = t(N)-t(1);

% analytical results

y = sym('y');

lambda_th = -40:1e-3:40;
Delta_f = @(x) sqrt(x.^2 + abs(A).^2);
Delta = Delta_f(lambda_th);
qc_theory = conj(A)/1i./lambda_th.*exp(-1i.*lambda_th.*t(N)) .*...
    (1-Delta./1i./lambda_th.*cot(Delta.*T)).^(-1);

a = @(x) (cos(Delta_f(x).*T) - 1i.*x./Delta_f(x).*sin(Delta_f(x).*T)).*exp(1i.*x.*T);
da(y) = diff(a(y),y);
b = @(x) -conj(A)./Delta_f(x) .* sin(Delta_f(x).*T).*exp(-1i.*x.*(t(N)+t(1)));

% find eig from analytical eq.
di = 1e-6;
imagAxis = [0:di:A+1].*1i;
equation = exp(2*1i*(T2-T1).*sqrt(imagAxis.^2+abs(A).^2)) - (imagAxis + sqrt( imagAxis.^2 + A^2 )) ./...
    (imagAxis - sqrt( imagAxis.^2 + A^2 ));
epsilon = di*10;
[~,ind] = find(real(equation)>-epsilon & real(equation)<epsilon & imag(equation)>-epsilon & imag(equation)<epsilon );
k=1; % compare position offset
j=1; % current position
i=1; % free position pointer;
while j+k<=length(ind)
    if ind(j)+k==ind(j+k)
        if j+k==length(ind)
            ind(i) = ind(j+ceil((k-1)/2));
        end
        k=k+1;
    else
        ind(i) = ind(j+ceil((k-1)/2));
        i=i+1;
        if j+k==length(ind)
            ind(i) = ind(end);
        end
        j=j+k;
        k=1;
    end
end
ind = ind(1:i);
lambdaj_th = imagAxis(ind(1:end-1)); % last one is the degenerate case
qd_th = b(lambdaj_th)./double(da(lambdaj_th));

% NFT numerical - Forward
[qc,lambda,qd,lambdaj] = NFT_search(q,t,'fw');

% Plots
figure(1)
plot(lambda_th,abs(qc_theory), lambda, abs(qc), '.')
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