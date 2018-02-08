%read time values

clear all %#ok<CLALL>
close all

load time_average_v2.mat

time_mean = squeeze(mean(time,1));
time_var = squeeze(var(time,1));

% compare
figure(1)
hold on
for i=1:length(discreteEigenvalues_TOT)
    plot(K,squeeze(time_mean(1,:,i)),'b-o',K,squeeze(time_mean(2,:,i)),'r-o')
end
hold off
grid on
title('Computational time')
xlabel('Number of samples of waveform')
ylabel('Computational time [s]')
legend('Algorithm 1','Algorithm 2','Location','Northwest')

figure(2)
hold on
for i=1:length(K)
    plot(1:length(discreteEigenvalues_TOT),squeeze(time_mean(1,i,:)),'b-o',1:length(discreteEigenvalues_TOT),squeeze(time_mean(2,i,:)),'r-o')
end
hold off
grid on
title('Computational time')
xlabel('Number of discrete eigenvalues')
ylabel('Computational time [s]')
legend('Algorithm 1','Algorithm 2','Location','Northwest')

figure(3)
hold on
for i=1:length(K)
%     err = squeeze(time_mean(1,i,:)./time_mean(2,i,:)).^2 .*...
%           squeeze(time_var(1,i,:)./time_mean(1,i,:).^2 +...
%                   time_var(2,i,:)./time_mean(2,i,:).^2 -...
%                   2.*);
    plot(1:length(discreteEigenvalues_TOT),squeeze(time_mean(1,i,:))./squeeze(time_mean(2,i,:)),'-o')
    %errorbar(1:length(discreteEigenvalues_TOT),squeeze(time_mean(1,i,:))./squeeze(time_mean(2,i,:)),...
    %    squeeze(time_var()))
    s{i} = sprintf('samples = %s^{%d}','2',log2(K(i)));
end
plot(0:1e-1:length(discreteEigenvalues_TOT),ones(1,length(0:1e-1:length(discreteEigenvalues_TOT))),'k--')
hold off
grid on
title('Computational time ratio Alg1 / Alg2')
xlabel('Number of discrete eigenvalues')
ylabel('Computational time ratio')
legend(s)