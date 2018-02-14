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
discreteEigenvalues = [0 + 0.5i];
discreteSpectrum = [-1i];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%

nPoints = 2^12; % minimum 2048
Rs = Fs/nPoints;

%INFT parameters
param.INFT.Tn        = 0.15/Rs;                              %[km]
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
param.NFT.tolAllEigenvaluesFound = 0.3; % energy threshold to fulfill for terminating the eigs search
param.NFT.computeDiscreteSpectrumEnabled = 1;
param.NFT.computeContinuousSpectrumEnabled = 1;
param.NFT.complexPlaneSearchArea = 3 * (max(0.1, max(real(discreteEigenvalues))) + ...
    1i*max(imag(discreteEigenvalues)));
param.NFT.mexEnabled = 1;
param.NFT.returnNFTParameterB = param.INFT.setNFTParameterB;

% Let's get the normalization parameters that we need to build the reference signal
inft = DiscreteINFT_v1(param.INFT);
inft.normalizationParameters(const.c/Fc);

% Compute waveform with Darboux transform
inft = DiscreteINFT_v1(param.INFT);
sigDarb = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);

nft_out = NFT_v8(param.NFT);

for noise_index = 1:numel(osnr_cycle)
    
    %OSNR parameters
    param.OSNR.OSNR = osnr_cycle(noise_index);               % [dB]
    osnr = OSNR_v1(param.OSNR);
    
    parfor n=1:realization % realizations
        
        % Add noise to the generated signal
        sigNoise = osnr.traverse(sigDarb);
        
        % Compute NFT until all eignevalues are found
        nft_out.traverse(sigNoise);
        E = nft_out.results.E;
        
        robolog('Discrete eigenvalues OUTPUT signal:');
        egDb(n,:) = nft_out.discreteEigenvalues(); %egDb
        robolog('Discrete spectrum');
        dsDb(n,:) = nft_out.discreteSpectrum(); %dsDb
        E_cont(n) = E.Ec;
        E_disc(n) = E.Ed;
        t = genTimeAxisSig(sigNoise,'central');
        E_sigNoise(n) = trapz(t,abs(get(sigNoise)).^2);
        tmp{noise_index,n} = sigNoise;
    end
    store_sigNoise{noise_index} = tmp{noise_index,1};
    % check how energy flows from discrete to continuous spectrum
    E_cont_tot(noise_index) = mean(E_cont);
    E_disc_tot(noise_index) = mean(E_disc);
    E_sig_tot(noise_index) = mean(E_sigNoise);
    
    figure(1)
    plot(real(egDb),imag(egDb),'b.',real(discreteEigenvalues),imag(discreteEigenvalues),'r.')
    grid on
    xlabel('Real')
    ylabel('Imag')
    legend('Noisy eigs','Reference')
    
    figure(2)
    plot(real(dsDb),imag(dsDb),'b.',real(discreteSpectrum),imag(discreteSpectrum),'r.')
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

figure;
show_osnr_level = 1;
t = genTimeAxisSig(store_sigNoise{show_osnr_level},'central');
subplot(2,2,[1,3])
plot(t,abs(get(store_sigNoise{show_osnr_level})).^2,t,abs(get(sigDarb)).^2)
xlabel('t')
ylabel('|x(t)|^2')
subplot(2,2,[2,4])
plot(t,angle(get(store_sigNoise{show_osnr_level}))./pi,t,angle(get(sigDarb))./pi)
xlabel('t')
ylabel('phase/\pi')

figure;
plot(osnr_cycle,E_cont_tot,'-o',osnr_cycle,E_disc_tot,'-o',...
    osnr_cycle,E_cont_tot+E_disc_tot,'k--o',osnr_cycle,E_sig_tot,'k-o');
xlabel('OSNR [dB]')
ylabel('Signal Energy [J]')
legend('Continuous spectrum','Discrete Spectrum',...
    'Total NFT enegry','Total signal energy','Location','Northwest')
%}