% Darboux ver1 VS ver2 testing.
% author: Nicola De Renzis
% date: 01/02/2018

clear all %#ok<CLALL>
close all

Fs = 1e12;
Fc = 193.4e12;

testType = 1;
% Test 1: Compare computational time. Generate a N-Solition from a given spectrum. A N-soliton has N eigenvalue on the imaginary axis.
% Test 2: Compare precision darboux VS darboux_simplified. Generate a N-Solition from a given spectrum. A N-soliton has N eigenvalue on the imaginary axis.
% Test 3: Compare results of both algorithms against theoretical results
%         with one eignevalue only

switch testType
    case 1
        %% Test 1: Generate a N-Solition from a given spectrum
        % A N-soliton has N eigenvalue on the imaginary axis.
        % compare computational time
        
        setpref('roboLog', 'logLevel', 1)
        
        method = {'darboux','darboux_simplified'};
        % Set the spectrum here
        discreteEigenvalues_TOT = [1:4:4*8] + 1i.*[1:6:6*8];
        discreteSpectrum_TOT = rand(1,numel(discreteEigenvalues_TOT))*500-25 ...
            + 1i.*(rand(1,numel(discreteEigenvalues_TOT))*500-30);
        
        %K = 2.^(10:2:18); % length of signals
        K = 30:5:100; % length of signals
        TIME_AVERAGE = 2000;
        %K = 2^10;
        for method_cycle = 1:length(method)
            for k=1:length(K)
                nPoints = K(k);
                Rs = Fs/nPoints;
                for n=1:length(discreteEigenvalues_TOT)
                    
                    fprintf('%d %d\n',length(K)-k,length(discreteEigenvalues_TOT)-n)
                    
                    discreteEigenvalues = discreteEigenvalues_TOT(1:n);
                    discreteSpectrum = discreteSpectrum_TOT(1:n);
                    
                    %INFT parameters
                    param.INFT.Tn        = 0.15/Rs;                              %[km]
                    param.INFT.gamma     = 1.27;                             %[/W/km]
                    param.INFT.D         = 17;                            %[ps/(nm km)]
                    param.INFT.method    = method{method_cycle};
                    param.INFT.Fc        = Fc;
                    param.INFT.nPoints   = nPoints;
                    param.INFT.setNFTParameterB = 0;
                    
                    % Let's get the normalization parameters that we need to build the reference signal
                    inft = DiscreteINFT_v1(param.INFT);
                    inft.normalizationParameters(const.c/Fc);
                    
                    % Compute waveform with Darboux transform
                    inft = DiscreteINFT_v1(param.INFT);
                    %pause(0.1); % relaxation time
                    tic
                    for i=1:TIME_AVERAGE
                        sigDarb = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);
                    end
                    tmp = toc;
                    time(method_cycle,k,n)=tmp/TIME_AVERAGE;
                    
                    clear param inft discreteEigenvalues discreteSpectrum
                end
            end
        end
        save time_average2.mat time
        
        %% compare
        close all
        figure(1)
        hold on
        for i=1:length(discreteEigenvalues_TOT)
            plot(K,squeeze(time(1,:,i)),'b-o',K,squeeze(time(2,:,i)),'r-o')
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
            plot(1:length(discreteEigenvalues_TOT),squeeze(time(1,i,:)),'b-o',1:length(discreteEigenvalues_TOT),squeeze(time(2,i,:)),'r-o')
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
            plot(1:length(discreteEigenvalues_TOT),squeeze(time(1,i,:))./squeeze(time(2,i,:)),'-o')
            %s{i} = sprintf('samples = %s^{%d}','2',log2(K(i)));
            s{i} = sprintf('%d',K(i));
        end
        plot(0:1e-1:length(discreteEigenvalues_TOT),ones(1,length(0:1e-1:length(discreteEigenvalues_TOT))),'k--')
        hold off
        grid on
        title('Computational time ratio Alg1 / Alg2')
        ylim([0 max(max(time(1,:,:)./time(2,:,:)))+0.5])
        xlabel('Number of discrete eigenvalues')
        ylabel('Computational time ratio')
        legend(s)
        
        
    case 2
        %% Test 2: Generate a N-Solition from a given spectrum
        % A N-soliton has N eigenvalue on the imaginary axis.
        % compare darboux VS darboux_simplified
        
        % Set the spectrum here
        discreteEigenvalues_TOT = [1:4:50*2] + 1i.*[1:6:50*3];
        discreteSpectrum_TOT = rand(1,50/2)*500-25 + 1i.*(rand(1,50/2)*500-30);
        
        K = 2.^(10:2:18); % length of signals
        error = zeros(length(K),length(discreteEigenvalues_TOT));
        for k=1:length(K)
            disp(K-k)
            nPoints = K(k);
            Rs = Fs/nPoints;
            for n=1:length(discreteEigenvalues_TOT)
                discreteEigenvalues = discreteEigenvalues_TOT(1:n);
                discreteSpectrum = discreteSpectrum_TOT(1:n);
                
                %INFT parameters
                param.INFT.Tn        = 0.15/Rs;                              %[km]
                param.INFT.gamma     = 1.27;                             %[/W/km]
                param.INFT.D         = 17;                            %[ps/(nm km)]
                param.INFT.method    = 'darboux';
                param.INFT.Fc        = Fc;
                param.INFT.nPoints   = nPoints;
                param.INFT.setNFTParameterB = 0;
                
                % ALGORITHM 1
                % Let's get the normalization parameters that we need to build the reference signal
                inft = DiscreteINFT_v1(param.INFT);
                inft.normalizationParameters(const.c/Fc);
                % Compute waveform with Darboux transform
                inft = DiscreteINFT_v1(param.INFT);
                sigDarb_v1 = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);
                
                clear inft
                
                % ALGORITHM 2
                param.INFT.method = 'darboux_simplified';
                % Let's get the normalization parameters that we need to build the reference signal
                inft = DiscreteINFT_v1(param.INFT);
                inft.normalizationParameters(const.c/Fc);
                % Compute waveform with Darboux transform
                inft = DiscreteINFT_v1(param.INFT);
                sigDarb_v2 = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);
                
                error(k,n) = mean(abs(get(sigDarb_v1) - get(sigDarb_v2)).^2); % difference between them
                
                clear param inft discreteEigenvalues discreteSpectrum
            end
        end
        
        %%
        close all
        
        figure(1)
        for i=1:length(K)
            semilogy(1:length(discreteEigenvalues_TOT),error(i,:))
            hold on
            s{i} = sprintf('samples = %s^{%d}','2',log2(K(i)));
        end
        hold off
        grid on
        ylim([1e-20 1e0])
        title('Error between Alg1 VS Alg2')
        xlabel('Number of discrete eigenvalues')
        ylabel('MSE')
        legend(s,'Location','Southeast')
        
        figure(2)
        for i=1:length(discreteEigenvalues_TOT)
            semilogy(K,error(:,i),'b')
            hold on
        end
        hold off
        grid on
        %ylim([1e-20 1e0])
        title('Error between Alg1 VS Alg2')
        xlabel('Number of samples in waveform')
        ylabel('MSE')
        
        
        
    case 3
        %% Test 3: Compare results of both algorithms against theoretical results
        % with one eignevalue only
        
        % Set the spectrum here
        discreteEigenvalues = 11+1*1i;
        discreteSpectrum = -1000+3000*1i;
        
        K = 2.^(10:2:18); % length of signals
        
        error_v1 = zeros(1,length(K));
        error_v2 = zeros(1,length(K));
        error_v12 = zeros(1,length(K));
        for k=1:length(K)
            nPoints = K(k);
            Rs = Fs/nPoints;
            
            %INFT parameters
            param.INFT.Tn        = 0.03/Rs;                              %[km]
            param.INFT.gamma     = 1.27;                             %[/W/km]
            param.INFT.D         = 17;                            %[ps/(nm km)]
            param.INFT.method    = 'darboux';
            param.INFT.Fc        = Fc;
            param.INFT.setNFTParameterB = 0;
            
            % THEORY
            inft = DiscreteINFT_v1(param.INFT);
            inft.normalizationParameters(const.c/Fc);
            A = imag(discreteEigenvalues) * 2;
            t = (-20:1/Fs/inft.Tn:20)*inft.Tn;
            %dt = t(2)-t(1);
            %t = t-dt/2;
            if param.INFT.setNFTParameterB == 0
                t0 = 1/A*log(abs(discreteSpectrum)/A);
                phase = 2*real(discreteEigenvalues).*t/inft.Tn+angle(discreteSpectrum) + pi/2; % Note that pi and -pi are the same phase (in figure)
            else
                t0 = 1/A*log(abs(discreteSpectrum));
                phase = 2*real(discreteEigenvalues).*t/inft.Tn+angle(discreteSpectrum) + pi;
            end
            x = A*sqrt(inft.Pn)*sech(A*(t/inft.Tn - t0)).*exp(-1i*phase);
            sig = signal_interface(x, struct('Rs', inf, 'Fs', Fs, 'Fc', Fc));
            param.INFT.nPoints = sig.L;
            lengths(k) = sig.L;
            
            % ALGORITHM 1
            % Let's get the normalization parameters that we need to build the reference signal
            inft = DiscreteINFT_v1(param.INFT);
            inft.normalizationParameters(const.c/Fc);
            % Compute waveform with Darboux transform
            inft = DiscreteINFT_v1(param.INFT);
            sigDarb_v1 = inft.traverse(discreteEigenvalues, discreteSpectrum, Fs/inft.nPoints);
            % ALGORITHM 2
            param.INFT.method = 'darboux_simplified';
            % Let's get the normalization parameters that we need to build the reference signal
            inft = DiscreteINFT_v1(param.INFT);
            inft.normalizationParameters(const.c/Fc);
            % Compute waveform with Darboux transform
            inft = DiscreteINFT_v1(param.INFT);
            sigDarb_v2 = inft.traverse(discreteEigenvalues, discreteSpectrum, Fs/inft.nPoints);
            
            error_v1(k) = mean(abs(get(sig) - get(sigDarb_v1)).^2);
            error_v2(k) = mean(abs(get(sig) - get(sigDarb_v2)).^2);
            error_v12(k) = mean(abs(get(sigDarb_v1) - get(sigDarb_v2)).^2);
            clear param inft
            
        end
        
        %% figures
        figure(1)
        loglog(lengths,error_v1,'-o',lengths,error_v2,'-o')
        legend('Ver 1','Ver 2')
        xlabel('nPoints')
        ylabel('MSE')
        grid on
        
        figure(2)
        loglog(lengths,error_v12,'-o')
        legend('Ver 1 VS Ver 2')
        xlabel('nPoints')
        ylabel('MSE')
        grid on
        
end % end switch
