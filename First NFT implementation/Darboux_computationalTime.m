% Darboux ver1 VS ver2 testing.
% author: Nicola De Renzis
% date: 01/02/2018

clear all %#ok<CLALL>
close all

Fs = 1e12;
Fc = 193.4e12;

testType = 4;
% Test 1: Compare computational time. Generate a N-Solition from a given spectrum. A N-soliton has N eigenvalue on the imaginary axis.
% Test 2: Compare precision darboux VS darboux_simplified. Generate a N-Solition from a given spectrum. A N-soliton has N eigenvalue on the imaginary axis.
% Test 3: Compare results of both algorithms against theoretical results
%         with one eignevalue only
% Test 4: Compare results of both algorithms against theoretical results
%         with N eignevalue only (may take long, scaling with number of
%         samples)

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
        discreteEigenvalues = 1+1*1i;
        discreteSpectrum = -1+3*1i;
        
        K = 2.^(10:2:18); % length of signals
        K = 5:5:30;
        error_v1 = zeros(1,length(K));
        error_v2 = zeros(1,length(K));
        error_v12 = zeros(1,length(K));
        for k=1:length(K)
            nPoints = K(k);
            Rs = Fs/nPoints;
            
            %INFT parameters
            param.INFT.Tn        = 0.15/Rs;                              %[km]
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
        
        
    case 4
        %         %% Test 4: Compare results of both algorithms against theoretical results
        %         % with N eignevalues (may take long)
        %         warning off
        %         setpref('roboLog', 'logLevel', 1)
        %
        %         N_average = 100;
        %         N_max_eigs = 6;
        %         samples = 40:1:150;
        %
        %         for average=1:N_average % average several simulations
        %             % Set the spectrum here
        %             discreteEigenvalues_TOT = rand(1,50/2)*5-5/2 + 1i.*(rand(1,50/2)*12+0.1);
        %             discreteSpectrum_TOT = rand(1,50/2)*500-25 + 1i.*(rand(1,50/2)*500-30);
        %             N_tot = min(N_max_eigs,numel(discreteEigenvalues_TOT));
        %
        %             ratio_Tn = 0.3;
        %             T1 = -20;
        %             T2 = -T1;
        %             samples = 40:1:200;
        %             K = (samples-1)./((T2-T1)*ratio_Tn);
        %
        %             error_v1 = zeros(N_tot,length(K));
        %             error_v2 = zeros(N_tot,length(K));
        %             error_v12 = zeros(N_tot,length(K));
        %             for n=1:N_tot
        %                 disp([N_average-average,N_tot-n])
        %                 discreteEigenvalues = discreteEigenvalues_TOT(1:n);
        %                 discreteSpectrum = discreteSpectrum_TOT(1:n);
        %
        %                 N = numel(discreteEigenvalues);
        %                 c = discreteSpectrum; % name only
        %                 zeta = discreteEigenvalues; % name only
        %
        %                 for i=1:length(K)
        %                     nPoints = K(i);
        %                     Rs = Fs/nPoints;
        %
        %                     %INFT parameters
        %                     param.INFT.Tn        = ratio_Tn/Rs;                  %[km]
        %                     param.INFT.gamma     = 1.27;                             %[/W/km]
        %                     param.INFT.D         = 17;                            %[ps/(nm km)]
        %                     param.INFT.method    = 'darboux';
        %                     param.INFT.Fc        = Fc;
        %                     param.INFT.setNFTParameterB = 0;
        %
        %                     % normalization of the resulting signal
        %                     Ld_SI = param.INFT.D * 1e3;
        %                     gamma_SI = param.INFT.gamma / 1e3;
        %                     D_SI = param.INFT.D / 1e6;
        %                     beta2_SI = (D_SI*(const.c/Fc)^2)/(2*pi*const.c);
        %                     Pn = abs(beta2_SI)/(gamma_SI*param.INFT.Tn^2);
        %
        %                     % THEORY
        %                     inft = DiscreteINFT_v1(param.INFT);
        %                     inft.normalizationParameters(const.c/Fc);
        %                     t = (T1:1/Fs/inft.Tn:T2)*inft.Tn;
        %                     dt = t(2)-t(1);
        %                     lambdas = @(r,z) sqrt(c(r)).*exp(1i.*zeta(r).*z); % index (j/k, time)
        %                     a = @(x,y,z) lambdas(x,z)*conj(lambdas(y,z))/(zeta(x)-conj(zeta(y))); % index (j,k,time)
        %                     b = @(x,y,z) conj(lambdas(x,z))*lambdas(y,z)/(conj(zeta(x))-zeta(y));
        %                     x = zeros(1,length(t));
        %                     for tau=1:length(t)
        %                         %disp(length(t)-tau);
        %                         t_curr = t(tau)/inft.Tn - dt/2/inft.Tn;
        %                         A_tilde = zeros(N,N);
        %                         B_tilde = zeros(N,N);
        %                         for j=1:N
        %                             %cj = c(j);
        %                             %zetaj = zeta(j);
        %                             for k=1:N
        %                                 A_tilde(j,k) = a(j,k,t_curr);
        %                                 B_tilde(j,k) = -b(j,k,t_curr);
        %                                 % A_tilde(j,k) = sqrt(cj)*conj(sqrt(c(k)))/...
        %                                 %               (zetaj-conj(zeta(k)))*...
        %                                 %               exp(1i*(zetaj-conj(zeta(k)))*t_curr);
        %                                 % B_tilde(j,k) = -sqrt(c(k))*conj(sqrt(cj))/...
        %                                 %               (-zeta(k)+conj(zetaj))*...
        %                                 %               exp(1i*(-zeta(k)+conj(zetaj))*t_curr);
        %                             end
        %                         end
        %                         A_tot = [eye(N),A_tilde;B_tilde,eye(N)];
        %                         b_tot = [zeros(1,N),conj(lambdas(1:N,t_curr))].';
        %                         psi = linsolve(A_tot,b_tot);
        %                         x(tau) = -2.*sum(b_tot(N+1:end).*psi(N+1:end))*1i;
        %                     end
        %                     % normalization of the resulting signal
        %                     x = x.*sqrt(Pn);
        %                     x(isnan(real(x))) = 1i.*imag(x(isnan(real(x))));
        %                     x(isnan(imag(x))) = real(x(isnan(imag(x))));
        %                     x(isnan(real(x))) = 1i.*imag(x(isnan(real(x))));
        %
        %                     sig = signal_interface(x, struct('Rs', inf, 'Fs', Fs, 'Fc', Fc));
        %                     param.INFT.nPoints = sig.L;
        %                     lengths(i) = sig.L;
        %
        %                     % ALGORITHM 1
        %                     % Let's get the normalization parameters that we need to build the reference signal
        %                     inft = DiscreteINFT_v1(param.INFT);
        %                     inft.normalizationParameters(const.c/Fc);
        %                     % Compute waveform with Darboux transform
        %                     inft = DiscreteINFT_v1(param.INFT);
        %                     sigDarb_v1 = inft.traverse(discreteEigenvalues, discreteSpectrum, Fs/inft.nPoints);
        %
        %                     % ALGORITHM 2
        %                     param.INFT.method = 'darboux_simplified';
        %                     % Let's get the normalization parameters that we need to build the reference signal
        %                     inft = DiscreteINFT_v1(param.INFT);
        %                     inft.normalizationParameters(const.c/Fc);
        %                     % Compute waveform with Darboux transform
        %                     inft = DiscreteINFT_v1(param.INFT);
        %                     sigDarb_v2 = inft.traverse(discreteEigenvalues, discreteSpectrum, Fs/inft.nPoints);
        %
        %                     error_v1(n,i) = mean(abs(get(sig) - get(sigDarb_v1)).^2);
        %                     error_v2(n,i) = mean(abs(get(sig) - get(sigDarb_v2)).^2);
        %                     error_v12(n,i) = mean(abs(get(sigDarb_v1) - get(sigDarb_v2)).^2);
        %                     %figure;plot(t,real(get(sigDarb_v1)),'g',t,imag(get(sigDarb_v1)),'g--',t,real(get(sigDarb_v2)),'b',t,imag(get(sigDarb_v2)),'b--',t,real(get(sig)),'r',t,imag(get(sig)),'r--')
        %                     clear param inft
        %
        %                 end
        %             end
        %             error_v1_average(average,:,:) = error_v1;
        %             error_v2_average(average,:,:) = error_v2;
        %             error_v12_average(average,:,:) = error_v12;
        %         end
        %         error_v1 = squeeze(mean(error_v1_average,1));
        %         error_v2 = squeeze(mean(error_v2_average,1));
        %         error_v12 = squeeze(mean(error_v12_average,1));
        %
        %         %% figures
        %         figure(1)
        %         cmap = colormap(parula(N_tot+2));
        %         for i=1:N_tot
        %             semilogy(lengths,error_v1(i,:),'Color',cmap(i,:),'Marker','o','Markersize',1)
        %             hold on
        %             semilogy(lengths,error_v2(i,:),'Color',cmap(i,:),'LineStyle','--','Marker','o','Markersize',1)
        %         end
        %         hold off
        %         title('Error from theoretical results')
        %         legend('Ver 1','Ver 2')
        %         xlabel('samples')
        %         ylabel('MSE')
        %         grid on
        %
        %         figure(2)
        %         semilogy(lengths,error_v12,'-o','Markersize',1)
        %         title('Ver 1 VS Ver 2 - Number of eigenvalues')
        %         legend(num2cell(num2str((1:N_tot).')))
        %         xlabel('samples')
        %         ylabel('MSE')
        %         grid on
        %
        %         figure(3)
        %         semilogy(lengths,error_v1./error_v2,'-o','Markersize',1)
        %         title('Ver 1 VS Ver 2 - Error ratio')
        %         legend(num2cell(num2str((1:N_tot).')))
        %         xlabel('samples')
        %         ylabel('MSE ratio')
        %         grid on
        %
        %         % reset robolog
        %         setpref('roboLog', 'logLevel', 5)
        
        
        %% Test 4: Compare results of both algorithms against theoretical results
        % with N eignevalues (may take long)
        warning off
        setpref('roboLog', 'logLevel', 1)
        
        N_average = 200;
        N_max_eigs = 6;
        
        Nss = 40:1:150;
        ratio_Tn = 0.05;
        T1 = -50*1e-12;
        T2 = -T1;
        Rs = 1/(T2-T1);
        
        for average=1:N_average % average several simulations
            % Set the spectrum here
            discreteEigenvalues_TOT = rand(1,50/2)*5-5/2 + 1i.*(rand(1,50/2)*12+0.1);
            discreteSpectrum_TOT = rand(1,50/2)*500-25 + 1i.*(rand(1,50/2)*500-30);
            %discreteEigenvalues_TOT = 1i*0.5;
            %discreteSpectrum_TOT = -1i;
            N_tot = min(N_max_eigs,numel(discreteEigenvalues_TOT));
            
            error_v1 = zeros(N_tot,length(Nss));
            error_v2 = zeros(N_tot,length(Nss));
            error_v12 = zeros(N_tot,length(Nss));
            error_v1_inter = zeros(N_tot,length(Nss));
            error_v2_inter = zeros(N_tot,length(Nss));
            for n=1:N_tot
                disp([N_average-average,N_tot-n])
                discreteEigenvalues = discreteEigenvalues_TOT(1:n);
                discreteSpectrum = discreteSpectrum_TOT(1:n);
                
                N = numel(discreteEigenvalues);
                c = discreteSpectrum; % name only
                zeta = discreteEigenvalues; % name only
                
                for i=1:length(Nss)
                    %Nss = 4001;
                    t = linspace(T1,T2,Nss(i));
                    %dt = t(2)-t(1);
                    %dt2 = 1/(Rs*(Nss(i)-1));
                    Fs = Rs*Nss(i);
                    
                    %INFT parameters
                    param.INFT.Tn        = ratio_Tn/Rs;                  %[km]
                    param.INFT.gamma     = 1.27;                             %[/W/km]
                    param.INFT.D         = 17;                            %[ps/(nm km)]
                    param.INFT.method    = 'darboux';
                    param.INFT.Fc        = Fc;
                    param.INFT.setNFTParameterB = 0;
                    param.INFT.nPoints = Nss(i);
                    
                    % normalization of the resulting signal
                    Ld_SI = param.INFT.D * 1e3;
                    gamma_SI = param.INFT.gamma / 1e3;
                    D_SI = param.INFT.D / 1e6;
                    beta2_SI = (D_SI*(const.c/Fc)^2)/(2*pi*const.c);
                    Pn = abs(beta2_SI)/(gamma_SI*param.INFT.Tn^2);
                    
                    % THEORY
                    inft = DiscreteINFT_v1(param.INFT);
                    inft.normalizationParameters(const.c/Fc);

                    lambdas = @(r,z) sqrt(c(r)).*exp(1i.*zeta(r).*z); % index (j/k, time)
                    a = @(x,y,z) lambdas(x,z).*conj(lambdas(y,z))./(zeta(x)-conj(zeta(y))); % index (j,k,time)
                    b = @(x,y,z) conj(lambdas(x,z)).*lambdas(y,z)./(conj(zeta(x))-zeta(y));
                    x = zeros(1,length(t));
                    for tau=1:length(t)
                        %disp(length(t)-tau);
                        t_curr = t(tau)/inft.Tn;% - dt/2/inft.Tn;
                        A_tilde = zeros(N,N);
                        B_tilde = zeros(N,N);
                        for j=1:N
                            %cj = c(j);
                            %zetaj = zeta(j);
                            A_tilde(j,:) = a(j,1:N,t_curr);
                            B_tilde(j,:) = -b(j,1:N,t_curr);
                            %{
                            %for k=1:N
                                %A_tilde(j,k) = a(j,k,t_curr);
                                %B_tilde(j,k) = -b(j,k,t_curr);
                                % A_tilde(j,k) = sqrt(cj)*conj(sqrt(c(k)))/...
                                %               (zetaj-conj(zeta(k)))*...
                                %               exp(1i*(zetaj-conj(zeta(k)))*t_curr);
                                % B_tilde(j,k) = -sqrt(c(k))*conj(sqrt(cj))/...
                                %               (-zeta(k)+conj(zetaj))*...
                                %               exp(1i*(-zeta(k)+conj(zetaj))*t_curr);
                            %end
                            %}
                        end
                        A_tot = [eye(N),A_tilde;B_tilde,eye(N)];
                        b_tot = [zeros(1,N),conj(lambdas(1:N,t_curr))].';
                        psi = linsolve(A_tot,b_tot);
                        x(tau) = -2.*sum(b_tot(N+1:end).*psi(N+1:end))*1i; % THIS i IS NOT SURE, BUT IT WORKS
                    end
                    % normalization of the resulting signal
                    x = x.*sqrt(Pn);
                    x(isnan(real(x))) = 1i.*imag(x(isnan(real(x))));
                    x(isnan(imag(x))) = real(x(isnan(imag(x))));
                    x(isnan(real(x))) = 1i.*imag(x(isnan(real(x))));
                    
                    sig = signal_interface(x, struct('Rs', Rs, 'Fs', Fs, 'Fc', Fc));
                    
                    % ALGORITHM 1
                    % Let's get the normalization parameters that we need to build the reference signal
                    inft = DiscreteINFT_v1(param.INFT);
                    inft.normalizationParameters(const.c/Fc);
                    % Compute waveform with Darboux transform
                    %inft = DiscreteINFT_v1(param.INFT);
                    sigDarb_v1 = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);
                    % interpolate_v1 (add a sample at the end since upsampled time is 2N-1)
                    x_tmp = interp1(t,get(sigDarb_v1).',linspace(T1,T2,Nss(i)*2-1));
                    x_int = [x_tmp(2:2:end),x_tmp(end)];
                    sigDarb_v1_inter = signal_interface(x_int, struct('Rs', Rs, 'Fs', Fs, 'Fc', Fc));
                    
                    % ALGORITHM 2
                    param.INFT.method = 'darboux_simplified';
                    % Let's get the normalization parameters that we need to build the reference signal
                    inft = DiscreteINFT_v1(param.INFT);
                    inft.normalizationParameters(const.c/Fc);
                    % Compute waveform with Darboux transform
                    %inft = DiscreteINFT_v1(param.INFT);
                    sigDarb_v2 = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);
                    % interpolate_v1 (see up for definition)
                    x_tmp = interp1(t,get(sigDarb_v2).',linspace(T1,T2,Nss(i)*2-1));
                    x_int = [x_tmp(2:2:end),x_tmp(end)];
                    sigDarb_v2_inter = signal_interface(x_int, struct('Rs', Rs, 'Fs', Fs, 'Fc', Fc));

                    error_v1(n,i) = mean(abs(get(sig) - get(sigDarb_v1)).^2);
                    error_v2(n,i) = mean(abs(get(sig) - get(sigDarb_v2)).^2);
                    error_v1_inter(n,i) = mean(abs(get(sig) - get(sigDarb_v1_inter)).^2);
                    error_v2_inter(n,i) = mean(abs(get(sig) - get(sigDarb_v2_inter)).^2);
                    error_v12(n,i) = mean(abs(get(sigDarb_v1) - get(sigDarb_v2)).^2);
                    %figure;plot(linspace(T1,T2,Nss(i)*2-1),real(x_tmp),'g--o',t,real(get(sigDarb_v2_inter)),'y-o',t,imag(get(sigDarb_v2_inter)),'y--o',t,real(get(sigDarb_v1_inter)),'r-o',t,imag(get(sigDarb_v1_inter)),'r--o',t,real(get(sigDarb_v2)),'b-o',t,imag(get(sigDarb_v2)),'b--o',t,real(get(sig)),'k-o',t,imag(get(sig)),'k--o');grid on
                    clear param inft
                    
                end
            end
            error_v1_average(average,:,:) = error_v1;
            error_v2_average(average,:,:) = error_v2;
            error_v1_average_inter(average,:,:) = error_v1_inter;
            error_v2_average_inter(average,:,:) = error_v2_inter;
            error_v12_average(average,:,:) = error_v12;
        end
        error_v1 = squeeze(mean(error_v1_average,1));
        error_v2 = squeeze(mean(error_v2_average,1));
        error_v12 = squeeze(mean(error_v12_average,1));
        error_v1_inter = squeeze(mean(error_v1_average_inter,1));
        error_v2_inter = squeeze(mean(error_v2_average_inter,1));
        
        %% figures
        figure(1)
        cmap = colormap(parula(N_tot+2));
        for i=1:N_tot
            semilogy(Nss,error_v1(i,:),'Color',cmap(i,:),'Marker','o','Markersize',1)
            hold on
            semilogy(Nss,error_v2(i,:),'Color',cmap(i,:),'LineStyle','--','Marker','o','Markersize',1)
        end
        hold off
        title('Error from theoretical results')
        legend('Ver 1','Ver 2')
        xlabel('samples')
        ylabel('MSE')
        grid on
        
        figure(2)
        cmap = colormap(parula(N_tot+2));
        for i=1:N_tot
            semilogy(Nss,error_v1_inter(i,:),'Color',cmap(i,:),'Marker','o','Markersize',1)
            hold on
            semilogy(Nss,error_v2_inter(i,:),'Color',cmap(i,:),'LineStyle','--','Marker','o','Markersize',1)
        end
        hold off
        title('Error from theoretical results - Interpolated')
        legend('Ver 1','Ver 2')
        xlabel('samples')
        ylabel('MSE')
        grid on
        
        figure(3)
        semilogy(Nss,error_v12,'-o','Markersize',1)
        title('Ver 1 VS Ver 2 - Number of eigenvalues')
        legend(num2cell(num2str((1:N_tot).')))
        xlabel('samples')
        ylabel('MSE')
        grid on
        
        figure(4)
        semilogy(Nss,error_v1./error_v2,'-o','Markersize',1)
        title('Ver 1 VS Ver 2 - Error ratio')
        legend(num2cell(num2str((1:N_tot).')))
        xlabel('samples')
        ylabel('MSE ratio')
        grid on
        
        % reset robolog
        setpref('roboLog', 'logLevel', 5)
        
end % end switch
