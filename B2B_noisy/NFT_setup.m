clearall
close all
clc
robolog off

Fs = 1e12;
Fc = 193.4e12;

osnr_cycle = 42:-4:14;
realization = 1e2;


%% Generate a waveform from a given spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
% Set the spectrum here
discreteEigenvalues = [0.5i]; % ordered by increasing both imaginary and real part (one after the other: [-1-1i, 1-1i, -1+1i, 1+1i])
discreteSpectrum = [-1i];
N = numel(discreteEigenvalues);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%

N_padding = 20;
discreteEigenvalues = [discreteEigenvalues,zeros(1,N_padding-N)];
discreteSpectrum = [discreteSpectrum,zeros(1,N_padding-N)];
egDb = zeros(realization,N_padding);
dsDb = zeros(realization,N_padding);

nPoints = 2^10; % minimum 2048smart
Rs = Fs/nPoints;

%INFT parameters
param.INFT.Tn        = 0.07/Rs;                              %[km]
param.INFT.gamma     = 1.27;                             %[/W/km]
param.INFT.D         = 17;                            %[ps/(nm km)]
param.INFT.method    = 'darboux_simplified';
param.INFT.Fc        = Fc;
param.INFT.nPoints   = nPoints;
param.INFT.setNFTParameterB = 0;

%NFT parameters
param.NFT.Tn        = param.INFT.Tn;                              %[km]
param.NFT.gamma     = param.INFT.gamma;                             %[/W/km]
param.NFT.D         = param.INFT.D;                               %[ps/(nm km)]
param.NFT.nPoints   = nPoints;
param.NFT.methodDiscreteSpectrum = 'TrapezFB';
param.NFT.tolUniqueEigenvalue = 0.05; % sensibility to locating the same eigenvalue
param.NFT.tolAllEigenvaluesFound = 0.1; % energy threshold to fulfill for terminating the eigs search
param.NFT.computeDiscreteSpectrumEnabled = 1;
param.NFT.computeContinuousSpectrumEnabled = 1;
param.NFT.complexPlaneSearchArea = 3 * (max(1, max(real(discreteEigenvalues(1:N)))) + ...
    1i*max(imag(discreteEigenvalues(1:N))));
param.NFT.mexEnabled = 1;
param.NFT.returnNFTParameterB = param.INFT.setNFTParameterB;

% Let's get the normalization parameters that we need to build the reference signal
%inft = DiscreteINFT_v1(param.INFT);
%inft.normalizationParameters(const.c/Fc);

% Compute waveform with Darboux transform
inft = DiscreteINFT_v1(param.INFT);
sigDarb = inft.traverse(discreteEigenvalues(1:N), discreteSpectrum(1:N), Rs);

nft_out = NFT_v8(param.NFT);
timeout = 10;

for noise_index = 1:numel(osnr_cycle)
    
    %OSNR parameters
    param.OSNR.OSNR = osnr_cycle(noise_index);               % [dB]
    osnr = OSNR_v1(param.OSNR);
    flag_start = 1;
    flag = 0;
    
    for n=1:realization % realizations
        
        % Add noise to the generated signal
        while flag_start || flag
            flag_start = 0;
            sigNoise = osnr.traverse(sigDarb);
            
            counter=0;
            toten=0;
            while (toten<0.8 || toten>1.001) && counter<timeout
                % Compute NFT until all eignevalues are found
                nft_out.traverse(sigNoise);
                E = nft_out.results.E;
                
                E_cont(n) = E.Ec;
                E_disc(n) = E.Ed;
                t = genTimeAxisSig(sigNoise,'central');
                E_sigNoise(n) = trapz(t./nft_out.Tn,abs(get(sigNoise)./sqrt(nft_out.Pn)).^2);
                tmp{noise_index,n} = sigNoise;
                
                toten = (E.Ec+E.Ed)./E_sigNoise(n);
                check=nft_out.discreteEigenvalues();
                if N>1 % only if there are more than 1 eigenvalue
                    if (~any(gradient(check)) || counter+1>timeout) && ~(toten<0.8 || toten>1.001) % check if toten condition is already not fulfilled
                        flag=1;
                        counter=1e16;
                    else
                        counter = counter+1;
                        flag=0;
                    end
                else
                    counter = counter+1;
                end
            end 
        end
        flag_start = 1;
        
        disp([numel(osnr_cycle)-noise_index,realization-n,counter-1])
        % order reults
        tmp_eigs = nft_out.discreteEigenvalues();
        tmp_amp = nft_out.discreteSpectrum();
        [~,index] = sort(real(tmp_eigs));
        tmp_eigs = tmp_eigs(index);
        tmp_amp = tmp_amp(index);
        [~,index] = sort(imag(tmp_eigs));
        tmp_eigs = tmp_eigs(index);
        tmp_amp = tmp_amp(index);
        
        % remember in which iteration a spurious eigenvalue was found
        if length(tmp_eigs)>N
            spurious{noise_index,n}.more =  length(tmp_eigs)-N;
            spurious{noise_index,n}.eigs = tmp_eigs;
            spurious{noise_index,n}.amp = tmp_amp;
        end
        
        egDb(n,:) = [tmp_eigs,zeros(1,N_padding-length(tmp_eigs))]; %egDb
        dsDb(n,:) = [tmp_amp,zeros(1,N_padding-length(tmp_amp))]; %dsDb
        
    end
    store_sigNoise{noise_index} = tmp{noise_index,1};
    % check how energy flows from discrete to continuous spectrum
    E_cont_tot(noise_index) = mean(E_cont);
    E_disc_tot(noise_index) = mean(E_disc);
    E_sig_tot(noise_index) = mean(E_sigNoise);
    
    figure(1)
    for i=1:numel(discreteEigenvalues)
        plot(real(egDb(:,i)),imag(egDb(:,i)),'.',real(discreteEigenvalues(i)),imag(discreteEigenvalues(i)),'ko')
        hold on
    end
    hold off
    grid on
    xlabel('Real')
    ylabel('Imag')
    legend('Noisy eigs','Reference')
    
    figure(2)
    for i=1:numel(discreteEigenvalues)
        plot(real(dsDb(:,i)),imag(dsDb(:,i)),'.',real(discreteSpectrum(i)),imag(discreteSpectrum(i)),'ko')
        hold on
    end
    hold off
    grid on
    xlabel('Real')
    ylabel('Imag')
    legend('Noisy amplitudes','Reference')
    
    %plotNFTConstellation_v1(egDb, dsDb , 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);
    mean_errorEigs(noise_index,:) = mean(egDb-discreteEigenvalues,1);
    mean_errorAmpl(noise_index,:) = mean(dsDb-discreteSpectrum,1);
    
    var_errorEigs(noise_index,:) = var(egDb,1);
    var_errorAmpl(noise_index,:) = var(dsDb,1);
    
end

%%{
%% Display the results
%plotNFTConstellation_v1(egDb, dsDb , 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);

figure(100);
show_osnr_level = 5;
t = genTimeAxisSig(store_sigNoise{show_osnr_level},'central');
subplot(2,2,[1,3])
plot(t/1e-12,abs(get(store_sigNoise{show_osnr_level})).^2,t/1e-12,abs(get(sigDarb)).^2)
grid on
xlabel('t [ps]')
ylabel('|x(t)|^2')
subplot(2,2,[2,4])
plot(t/1e-12,angle(get(store_sigNoise{show_osnr_level}))./pi,t/1e-12,angle(get(sigDarb))./pi)
grid on
xlabel('t [ps]')
ylabel('phase/\pi')

figure(101);
plot(osnr_cycle,E_cont_tot,'-o',osnr_cycle,E_disc_tot,'-o',...
    osnr_cycle,E_cont_tot+E_disc_tot,'k--o',osnr_cycle,E_sig_tot,'k-x');
xlabel('OSNR [dB]')
ylabel('Signal Energy [J]')
legend('Continuous spectrum','Discrete Spectrum',...
    'Total NFT enegry','Total signal energy','Location','Northeast')
title('Energy')
grid on

figure(102)
for i=1:numel(discreteEigenvalues)
    errorbar(osnr_cycle,abs(mean_errorEigs(:,i)),sqrt(abs(var_errorEigs(:,i))))
    hold on
end
hold off
title(['Averaged over ',num2str(realization),' iterations'] )
ylabel(' $| E [ \lambda - \hat{\lambda} ] |$','Interpreter','Latex')
xlabel('OSNR [dB]')
grid on

figure(103)
for i=1:numel(discreteEigenvalues)
    errorbar(osnr_cycle,abs(mean_errorAmpl(:,i)),sqrt(abs(var_errorAmpl(:,i))))
    hold on
end
hold off
title(['Averaged over ',num2str(realization),' iterations'] )
ylabel(' $| E [ Q_d(\lambda) - \hat{Q}_d(\hat{\lambda}) ] |$','Interpreter','Latex')
xlabel('OSNR [dB]')
grid on

%}