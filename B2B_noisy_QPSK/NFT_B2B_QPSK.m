clearall
clc
close all
robolog off

Fs = 1e12;
Fc = 193.4e12;

osnr_cycle = 18:-4:6;
%osnr_cycle = 10;
realization = 1*1e3;

trainModel = 0;

if trainModel
    %%%%%%%%%%%%%%%%%%%%%%% TRAIN AND LABELING SESSION %%%%%%%%%%%%%%%%%%%%%%%%
    disp('labeling...')
    discreteEigenvalues = [0.5*1i 1.5*1i]; % ordered by increasing both imaginary and real part (one after the other: [-1-1i, 1-1i, -1+1i, 1+1i])
    trainSpectrum = [-100i +100i];
    osnrTrain = 8;
    model = createModel(discreteEigenvalues,trainSpectrum,osnrTrain,1e3);
else
    load model_SVM_Linear.mat
end

%% Generate a waveform from a given spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
% Set the spectrum here
discreteEigenvalues = [0.5*1i 1.5*1i]; % ordered by increasing both imaginary and real part (one after the other: [-1-1i, 1-1i, -1+1i, 1+1i])
% Buelow way
multiplier = 2; % 2 -> 2x BPSK pi/2 shifted
                % 1 -> 2x QPSK pi/4 shifted
discreteSpectrum = [exp(1i*pi*(0:0.5*multiplier:1.5)); exp(1i*pi*(0:0.5*multiplier:1.5))];
discreteSpectrum(1,:) = discreteSpectrum(1,:)*exp(1i*0.25*multiplier*pi); % all possible compibations
discreteSpectrum = discreteSpectrum.';

[D,N] = size(discreteSpectrum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%

N_padding = 10;
discreteEigenvalues = [discreteEigenvalues,zeros(1,N_padding-N)];
%discreteSpectrum = [discreteSpectrum.',zeros(D,N_padding-N)];
egDb = zeros(realization,N);
dsDb = zeros(realization,N);

nPoints = 2^10; % minimum 2048
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


spurious_counter = 1;
flag_spurious = 0;
nft_out = NFT_v8(param.NFT);
timeout = 10;
totenMin = 0.8;
%%%%%%%%%%% INIT EMPTY VECTORS %%%%%%%%%%%%%
constSent = zeros(realization,D);
currentDiscreteSpectrum = zeros(realization,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for noise_index = 1:numel(osnr_cycle)
    
    %OSNR parameters
    param.OSNR.OSNR = osnr_cycle(noise_index);               % [dB]
    osnr = OSNR_v1(param.OSNR);
    flag_start = 1;
    flag = 0;
    counter_onlyOneFound = 0;
    iter_onylOne = [];
    
    for n=1:realization % realizations
        
        % choose random values to send from the allowed constellation
        constSent(n,:) = [randi(D),randi(D)];
        currentConst = constSent(n,:);
        for column = 1:N
            col = discreteSpectrum(:,column);
            currentDiscreteSpectrum(n,column) = col(currentConst(column));
        end
        
        % Compute waveform with Darboux transform
        inft = DiscreteINFT_v1(param.INFT);
        sigDarb = inft.traverse(discreteEigenvalues(1:N), currentDiscreteSpectrum(n,:), Rs);
        
        % Add noise to the generated signal
        while flag_start || flag
            flag_start = 0;
            sigNoise = osnr.traverse(sigDarb);
            
            counter=0;
            toten=0;
            while (toten<totenMin || toten>1.001) && counter<timeout
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
                    if (~any(gradient(check)) || counter+1>timeout) && ~(toten<totenMin || toten>1.001) % check if toten condition is already not fulfilled
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
        
        % order reults and apply model
        tmp_eigs = nft_out.discreteEigenvalues();
        tmp_amp = nft_out.discreteSpectrum();
        
        if numel(tmp_eigs)<N
           tmp_eigs = [tmp_eigs,tmp_eigs];
           tmp_amp = [tmp_amp,tmp_amp];
           counter_onlyOneFound = counter_onlyOneFound+1;
           iter_onylOne(counter_onlyOneFound+1) = n;
        end
        
        [label_tmp,loglike] = model_SVM_Linear.predictFcn([real(tmp_eigs);imag(tmp_eigs)].');
        if numel(label_tmp)>N
            loglike(loglike<0) = Inf;
            [~,minLoc] = min(loglike);
            label(n,:) = label_tmp(minLoc);
            tmp_eigs = tmp_eigs(minLoc);
            tmp_amp = tmp_amp(minLoc);
        else
            label(n,:) = label_tmp;
        end
        
        %{
        % remember in which iteration a spurious eigenvalue was found
        if length(tmp_eigs)>N
            spurious{spurious_counter}.more = length(tmp_eigs)-N;
            spurious{spurious_counter}.eigs = tmp_eigs;
            spurious{spurious_counter}.amp = tmp_amp;
            spurious{spurious_counter}.E = E;
            spurious{spurious_counter}.OSNR = osnr_cycle(noise_index);
            spurious_counter = spurious_counter+1;
            flag_spurious = 1;
        end
        %}
        egDb(n,:) = tmp_eigs; %egDb
        dsDb(n,:) = tmp_amp(label(n,:)); %dsDb
        
    end
    store_sigNoise{noise_index} = tmp{noise_index,1};
    % check how energy flows from discrete to continuous spectrum
    E_cont_tot(noise_index) = mean(E_cont);
    E_disc_tot(noise_index) = mean(E_disc);
    E_sig_tot(noise_index) = mean(E_sigNoise);
    
    % some figures
    figure(1)
    for i=1:N
        plot(real(egDb(label==i)),imag(egDb(label==i)),'.',real(discreteEigenvalues(i)),imag(discreteEigenvalues(i)),'ko')
        hold on
    end
    hold off
    grid on
    xlabel('Real')
    ylabel('Imag')
    legend('Noisy eigs','Reference')
    
    figure(2)
    for i=1:N
        plot(real(dsDb(:,i)),imag(dsDb(:,i)),'.',real(discreteSpectrum(:,min(N,i))),imag(discreteSpectrum(:,min(N,i))),'ko')
        hold on
    end
    hold off
    grid on
    xlabel('Real')
    ylabel('Imag')
    legend('Noisy amplitudes','Reference')
    
    %plotNFTConstellation_v1(egDb, dsDb , 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);
    %mean_errorEigs(noise_index,:) = mean(egDb-discreteEigenvalues,1);
    %mean_errorAmpl(noise_index,:) = mean(dsDb-discreteSpectrum,1);
    
    %var_errorEigs(noise_index,:) = var(egDb,1);
    %var_errorAmpl(noise_index,:) = var(dsDb,1);
    
    store{noise_index}.eigs = egDb;
    store{noise_index}.ampl = dsDb;
    store{noise_index}.osnr = osnr_cycle(noise_index);
    store{noise_index}.eigs_original = discreteEigenvalues(1:N);
    store{noise_index}.ampl_original = discreteSpectrum;
    if flag_spurious
        store{noise_index}.spurious = spurious;
    end
    
    %decode and check
    positions = 1:realization;
    for const=1:N % constellation are as many as the eigenvalues are
        recAmpl = dsDb(:,const);
        if const == 2
            recAmpl = recAmpl.*exp(1i*multiplier*pi/4); % rotate the ones on the axis
        end
        switch multiplier
            case 1
                constRec(positions(real(recAmpl)>0 & imag(recAmpl)>0),const) = 1;
                constRec(positions(real(recAmpl)<=0 & imag(recAmpl)>0),const) = 2;
                constRec(positions(real(recAmpl)<=0 & imag(recAmpl)<=0),const) = 3;
                constRec(positions(real(recAmpl)>0 & imag(recAmpl)<=0),const) = 4;
            case 2
                constRec(positions(imag(recAmpl)>0),const) = 1;
                constRec(positions(imag(recAmpl)<=0),const) = 2;
        end
    end
    error_mat = constRec - constSent;
    SER(noise_index) = sum(length(find(error_mat)))/N/realization;
    store{noise_index}.error_mat = constRec - constSent;
    store{noise_index}.SER =  SER(noise_index);
    oneFoundCounter(noise_index) = counter_onlyOneFound;
    SER_woOneEigFound(noise_index) =...
        (sum(length(find(error_mat(setdiff(positions,iter_onylOne))))))/...
        (N*(realization-counter_onlyOneFound));
end

%%{
%% Display the results
%plotNFTConstellation_v1(egDb, dsDb , 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);

figure(99)
semilogy(osnr_cycle,SER)

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
%{
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

name = sprintf('BPSK_%drealiz_%2.1fi.mat',realization,imag(discreteEigenvalues(N)));
save(name, 'store')
clear name