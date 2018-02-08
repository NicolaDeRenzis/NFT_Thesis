clearall
close all
clc

Fs = 1e12;
Fc = 193.4e12;

testType = 2;
caseType = 2;

% Test 1: Generate one soliton in time with arbitrary phase, compute NFT and then INFT and compare
    % Case 1: 1-soliton
    % Case 2: 2-soliton

% Test 2: Generate a 1-Solition from a given spectrum and compute NFT to compare with source spectrum
    % Case 1: use 'darboux' as INFT algorithm
    % Case 2: use 'darboux_v2' as INFT algorithm
% Test 3: Generate a 2-Solition from a given spectrum, compute NFT and compare spectrum
    % Case 1: use 'darboux' as INFT algorithm
    % Case 2: use 'darboux_v2' as INFT algorithm
% Test 4: Generate a 1-Solition for a full 2x QPSK constellation (compare time signal)
% Test 5: Generate a 1-Solition for a full 2x QPSK constellation (compare spectrum)
% Test 6: Generate a 2-Solitions for a full 2x QPSK constellation
switch testType
    case 1
        switch caseType
            %Figures:
            %1 Generated and darboux waveform
            %----
            %2 Input vs Darboux signal NFT eigenvalues
            %2 Input vs Darboux signal NFT spectrum
            
            %% Test:
            % Generate one soliton, compute NFT and then INFT and compare the two waveforms
            % in order to verify that the INFT works correctly
            
            % If A = B is a fundamental 1-soliton. Increasing A we increase the maximum amplitude and so the
            % power, the width of the sech changes as well.
            % If we increase the amplitude without changing the width, the sech will not be a fundamental
            % soliton any more for the given channel.
            % We can also set the phase
            case 1
                % CaseType 1 - 1-soliton (fundamental)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
                % Generate one 1-soliton
                A = 1;
                B = A;
                phase = 0.25*pi;
                t0u = 100*1/Fs;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%
            case 2
                % CaseType 2 - 2-soliton
                % RELATION BETWEEN A, B and eigenvalue
                % B = imag(discreteEigenvalues(2) - discreteEigenvalues(1));
                % A = imag(max(discreteEigenvalues)) + B/2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
                % Generate one 1-soliton
                % Works if abs(discreteSpectrum) < 1
                A = 1.2;
                B = 0.6;
                phase = 0.5*pi;
                t0u = 120*1/Fs;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%
            case 3
                % CaseType 3 - 3-soliton
                % RELATION BETWEEN A, B and eigenvalue
                % B = imag(discreteEigenvalues(2) - discreteEigenvalues(1));
                % A = imag(max(discreteEigenvalues)) + B/2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
                % Generate one 1-soliton
                % Works if abs(discreteSpectrum) < 1
                A = 1.8;
                B = 0.6;
                phase = 0.5*pi;
                t0u = 120*1/Fs;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%
        end
        
        %INFT parameters
        
        param.INFT.Tn        = 40*1/Fs;                              %[km]
        param.INFT.gamma     = 1.27;                             %[/W/km]
        param.INFT.D         = 17;                            %[ps/(nm km)]
        param.INFT.method    = 'darboux';
        param.INFT.Fc        = Fc;
        %NFT parameters
        param.NFT.Tn        = param.INFT.Tn;                              %[km]
        param.NFT.gamma     = param.INFT.gamma;                             %[/W/km]
        param.NFT.D         = param.INFT.D;                               %[ps/(nm km)]
        param.NFT.nPoints   = 2^7;
        param.NFT.method    = 'LP';                          % FCD ABL N-ABL LP
        param.NFT.methodDiscreteSpectrum = 'ABL';
        param.NFT.computeDiscreteSpectrumEnabled = 1;
        param.NFT.mexEnabled = 1;
        
        % Enable this to see the 3d plot of the 'a' coefficient an visually check that the NFT spectrum is
        % correct
        param.NFT.computeNftCoefficientsOnComplexPlaneEnabled = 0;
        
        % Let's get the NFT normalization parameters, we need them to build the proper soliton for the channel
        nft_in = NFT_v8(param.NFT);
        nft_in.normalizationParameters(const.c/Fc);
        
        t = (-12:1/Fs/nft_in.Tn:12)*nft_in.Tn;
        x = A*sqrt(nft_in.Pn)*sech(B*(t - t0u)/nft_in.Tn)*exp(1i*phase);
        sig = signal_interface(x, struct('Rs', inf, 'Fs', Fs, 'Fc', Fc));
        
        % Limit the search area of the eigenvalue, since we know where it is supposed to be.
        % The top right corner of the search rectangle should be specified.
        param.NFT.complexPlaneSearchArea = 0.05 + 1i*A; % This optimize the searching of the eigenvalue
        
        % Compute NFT and obtain its spectrum to feed it to the INFT
        nft_in = NFT_v8(param.NFT);
        nft_in.traverse(sig);
        discreteEigenvalues = nft_in.discreteEigenvalues();
        discreteSpectrum = nft_in.discreteSpectrum();
        
        % Now we compute the INFT with Darboux transform
        param.INFT.nPoints = sig.L;
        inft = DiscreteINFT_v1(param.INFT);
        inft.traverse(discreteEigenvalues, discreteSpectrum, Fs/param.INFT.nPoints);
        sigDarb = inft.timeDomainSignal();
        
        %% Display results
        % Plot input and output waveform
        pabsang(sig, sigDarb);
        legend('Source signal', 'Reconstructed signal');
        
        %% Extra plot and things
        disp('Press a button to show extra information...');
        pause;
        
        % Plot the 3d plot
        if param.NFT.computeNftCoefficientsOnComplexPlaneEnabled
            nft_in.plotNFTCoefficients3D;
        end
        
        % Compute NFT of waveform obtained from Darboux transform
        nft_out = NFT_v8(param.NFT);
        nft_out.traverse(sigDarb);
        
        try
            % Display textual results
            format compact;
            robolog('Discrete eigenvalues INPUT signal:');
            egRef = nft_in.discreteEigenvalues(); egRef
            robolog('Discrete spectrum');
            dsRef = nft_in.discreteSpectrum(); dsRef
        end
        
        try
            robolog('Discrete eigenvalues OUTPUT signal:');
            eg = nft_out.discreteEigenvalues(); eg
            robolog('Discrete spectrum');
            ds = nft_out.discreteSpectrum(); ds
        end
        
        try
            plotNFTConstellation_v1(eg, ds, 'refEigenvalues', egRef, 'refSpectrum', dsRef);
        end
        
    case 2
        %% Test 2: Generate a 1-Solition from a given spectrum
        % We compute NFT(INFT(eig)) and compare the constellations
        
        % We also generate the soliton in time manually and compare the
        % waveforms
        
        %Figures:
        %1 Input vs Recovered eigenvalues
        %2 Input vs Recovered discrete spectral amplitude
        %----
        %3 Eigenvalues recovered from manually generated signal
        %4 Spectral amplitudes recovered from manually generated signal
        %5 Darboux and manual signal waveforms
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
        % Set the spectrum here
        discreteEigenvalues = 3+1i*2;
        discreteSpectrum = 10 - 2*1i;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%
        
        %INFT parameters
        param.INFT.Tn        = 40*1/Fs;                              %[km]
        param.INFT.gamma     = 1.27;                             %[/W/km]
        param.INFT.D         = 17;                            %[ps/(nm km)]
        switch caseType
            case 1
                param.INFT.method    = 'darboux';       % possible values {'darboux','darboux_v2'}
            case 2
                param.INFT.method    = 'darboux_v2';    % possible values {'darboux','darboux_v2'}
        end
        param.INFT.nPoints   = 0; % Dummy value, real value is assigned below
        param.INFT.Fc        = Fc;
        param.INFT.setNFTParameterB = 0;
        
        %NFT parameters
        param.NFT.Tn        = param.INFT.Tn;                              %[km]
        param.NFT.gamma     = param.INFT.gamma;                             %[/W/km]
        param.NFT.D         = param.INFT.D;                               %[ps/(nm km)]
        param.NFT.nPoints   = 2^10;
        param.NFT.method    = 'TrapezFB';                          % FCD ABL N-ABL LP
        param.NFT.methodDiscreteSpectrum = 'TrapezFB';
        param.NFT.computeDiscreteSpectrumEnabled = 1;
        param.NFT.complexPlaneSearchArea = [1.1.*real(discreteEigenvalues) + 1i*1.1*imag(discreteEigenvalues)];
        param.NFT.mexEnabled = 1;
        
        % Let's get the normalization parameters that we need to build the reference signal
        inft = DiscreteINFT_v1(param.INFT);
        inft.normalizationParameters(const.c/Fc);
        
        % Autodetect A and B from the spectrum, theoretically we should be able to detect the phase as well
        %1-soliton
        A = imag(discreteEigenvalues) * 2;
        B = A;
        
        % Generate an 1-soliton manually (not needed in this case, we only use it to use the same number of
        % points in the Darboux signal)
        t = (-20:1/Fs/inft.Tn:20)*inft.Tn;
        if param.INFT.setNFTParameterB == 0
            t0 = 1/A*log(abs(discreteSpectrum)/A);
            phase = 2*real(discreteEigenvalues).*t/inft.Tn+angle(discreteSpectrum) + pi/2; % Note that pi and -pi are the same phase (in figure)
        else
            t0 = 1/A*log(abs(discreteSpectrum));
            phase = 2*real(discreteEigenvalues).*t/inft.Tn+angle(discreteSpectrum) + pi;
        end
        x = A*sqrt(inft.Pn)*sech(B*(t/inft.Tn - t0)).*exp(-1i*phase);
        sig = signal_interface(x, struct('Rs', inf, 'Fs', Fs, 'Fc', Fc));
        
        % Compute waveform with Darboux transform
        param.INFT.nPoints = sig.L;
        inft = DiscreteINFT_v1(param.INFT);
        sigDarb = inft.traverse(discreteEigenvalues, discreteSpectrum, Fs/inft.nPoints);
        
        nft_out = NFT_v8(param.NFT);
        nft_out.traverse(sigDarb);
        
        robolog('Discrete eigenvalues OUTPUT signal:');
        egDb = nft_out.discreteEigenvalues(); egDb
        robolog('Discrete spectrum');
        dsDb = nft_out.discreteSpectrum(); dsDb
        
        %% Display the results
        plotNFTConstellation_v1(egDb, dsDb, 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);
        
        %% Extra plot and things
        disp('Press a button to show extra information...');
        pause;
        
        % Compute NFT of manually generated soliton. It must have the same eigenvalues as the one set, or the
        % formula to generate it is wrong.
        nft_in = NFT_v8(param.NFT);
        nft_in.traverse(sig);
        
        robolog('Discrete eigenvalues INPUT signal:');
        eg = nft_in.discreteEigenvalues(); eg
        robolog('Discrete spectrum');
        ds = nft_in.discreteSpectrum(); ds
        plotNFTConstellation_v1(eg, ds, 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);
        
        % Plot the results
        figure;
        inft.timeDomainSignal();
        subplot(1,2,1); hold on;
        plot(t, abs(x), 'r--');
        subplot(1,2,2); hold on;
        plot(t, angle(x)/pi, 'r--');
    case 3
        %% Test 3: Generate a 2-Solition from a given spectrum
        % A 2-soliton has two eigenvalue on the imaginary axis. The waveform in time is a sech with A != B
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% TUNE PARAM %%%%%%%%%%%%%%%%%%%%%
        % Set the spectrum here
        discreteEigenvalues = [0 + 1.5i 0 + 0.5i];
        discreteSpectrum = [5*1i, 5-1i*15];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% END PARAM %%%%%%%%%%%%%%%%%%%%%
        
        nPoints = 2048; % Reduce this and you'll soon lose precision!
        %nPoints = 2^12;
        Rs = Fs/nPoints;
        
        %INFT parameters
        param.INFT.Tn        = 0.03/Rs;                              %[km]
        param.INFT.gamma     = 1.27;                             %[/W/km]
        param.INFT.D         = 17;                            %[ps/(nm km)]        switch caseType
        switch caseType
            case 1
                param.INFT.method    = 'darboux';       % possible values {'darboux','darboux_v2'}
            case 2
                param.INFT.method    = 'darboux_v2';    % possible values {'darboux','darboux_v2'}
        end
        param.INFT.Fc        = Fc;
        param.INFT.nPoints   = nPoints;
        param.INFT.setNFTParameterB = 0;
        
        %NFT parameters
        param.NFT.Tn        = param.INFT.Tn;                              %[km]
        param.NFT.gamma     = param.INFT.gamma;                             %[/W/km]
        param.NFT.D         = param.INFT.D;                               %[ps/(nm km)]
        param.NFT.nPoints   = 2^10;
        param.NFT.method    = 'TrapezFB';                          % FCD ABL N-ABL LP
        param.NFT.methodDiscreteSpectrum = 'TrapezFB';
        param.NFT.computeDiscreteSpectrumEnabled = 1;
        param.NFT.complexPlaneSearchArea = 1.1 * (max(0.1, max(real(discreteEigenvalues))) + ...
            1i*max(imag(discreteEigenvalues)));
        param.NFT.mexEnabled = 1;
        param.NFT.returnNFTParameterB = param.INFT.setNFTParameterB;
        % Let's get the normalization parameters that we need to build the reference signal
        inft = DiscreteINFT_v1(param.INFT);
        inft.normalizationParameters(const.c/Fc);
        
        % Compute waveform with Darboux transform
        inft = DiscreteINFT_v1(param.INFT);
        sigDarb = inft.traverse(discreteEigenvalues, discreteSpectrum, Rs);
        
        % Compute NFT until all eignevalues are found. A little random.
        toten = 0;
        while toten < 0.9
            nft_out = NFT_v8(param.NFT);
            nft_out.traverse(sigDarb);
            E = nft_out.results.E;
            toten = (E.Ec + E.Ed)/E.ETot;
        end
        
        robolog('Discrete eigenvalues OUTPUT signal:');
        egDb = nft_out.discreteEigenvalues(); egDb
        robolog('Discrete spectrum');
        dsDb = nft_out.discreteSpectrum(); dsDb
        
        %% Display the results
        plotNFTConstellation_v1(egDb, dsDb , 'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum);
        
        %% Extra plot and things
        disp('Press a button to show extra information...');
        pause;
        
        % Plot the results
        figure;
        inft.timeDomainSignal();
        
    case 4
        %% Test 4: Generate a 1-Solition for a full 2x QPSK constellation (compare time signal)
        
        %INFT parameters
        param.INFT.Tn        = 40*1/Fs;                              %[km]
        param.INFT.gamma     = 1.27;                             %[/W/km]
        param.INFT.D         = 17;                            %[ps/(nm km)]
        param.INFT.method    = 'darboux';
        param.INFT.Fc        = Fc;
        
        % Set the spectrum here
        discreteEigenvalues = 1i*[0.3 0.6].';
        %discreteSpectrum = -2*(discreteEigenvalues)*exp(1i*pi*[0:0.5:1.5]);
        discreteSpectrum = [exp(1i*pi*[0:0.5:1.5]); exp(1i*pi*[0:0.5:1.5])];
        discreteSpectrum(1,:) = discreteSpectrum(1,:)*exp(1i*0.25*pi);
        
        for i=1:length(discreteEigenvalues)
            for j=1:length(discreteSpectrum)
                % Let's get the normalization parameters that we need to build the reference signal
                inft = INFT_v2(param.INFT);
                inft.normalizationParameters(const.c/Fc);
                
                % Autodetect A and B from the spectrum, theoretically we should be able to detect the phase as well
                %1-soliton
                A = imag(discreteEigenvalues(i)) * 2;
                B = A;
                
                % Generate an 1-soliton manually
                phase = angle(discreteSpectrum(i,j)) + pi/2;
                t0 = 1/A*log(abs(discreteSpectrum(i,j))/A);
                t = (-12:1/Fs/inft.Tn:12)*inft.Tn;
                x = A*sqrt(inft.Pn)*sech(B*(t/inft.Tn - t0))*exp(-1i*phase);
                sig = signal_interface(x, struct('Rs', inf, 'Fs', Fs, 'Fc', Fc));
                
                % Compute waveform with Darboux transform
                inft = DiscreteINFT_v1(param.INFT);
                sigDarb = inft.traverse(discreteEigenvalues(i), discreteSpectrum(i,j), Fs/inft.nPoints);
                
                % Plot the results
                plotNFTConstellation_v1(discreteEigenvalues(i), discreteSpectrum(i,j));
                figure;
                inft.timeDomainSignal();
                subplot(1,2,1); hold on;
                plot(t, abs(x), 'r--');
                subplot(1,2,2); hold on;
                plot(t, angle(x)/pi, 'r--');
            end
        end
    case 5
        %% Test 5: Generate a 1-Solition for a full 2x QPSK constellation (compare spectra)
        
        nPoints = 2048; % Reduce this and you'll soon lose precision!
        Rs = Fs/nPoints;
        
        %INFT parameters
        param.INFT.Tn        = 0.03/Rs;                              %[km]
        param.INFT.gamma     = 1.27;                             %[/W/km]
        param.INFT.D         = 17;                            %[ps/(nm km)]
        param.INFT.method    = 'darboux';
        param.INFT.Fc        = Fc;
        param.INFT.nPoints   = nPoints;
        
        % Set the spectrum here
        discreteEigenvalues = 1i*[0.3 0.6].';
        %discreteSpectrum = -2*(discreteEigenvalues)*exp(1i*pi*[0:0.5:1.5]);
        discreteSpectrum = [exp(1i*pi*[0:0.5:1.5]); exp(1i*pi*[0:0.5:1.5])];
        discreteSpectrum(1,:) = discreteSpectrum(1,:)*exp(1i*0.25*pi);
        
        %NFT parameters
        param.NFT.Tn        = param.INFT.Tn;                              %[km]
        param.NFT.gamma     = param.INFT.gamma;                             %[/W/km]
        param.NFT.D         = param.INFT.D;                               %[ps/(nm km)]
        param.NFT.nPoints   = 2^10;
        param.NFT.method    = 'ABL';                          % FCD ABL N-ABL LP
        param.NFT.methodDiscreteSpectrum = 'ABL';
        param.NFT.computeDiscreteSpectrumEnabled = 1;
        param.NFT.complexPlaneSearchArea = 1.1 * (max(0.1, max(real(discreteEigenvalues))) + ...
            1i*max(imag(discreteEigenvalues)));
        param.NFT.mexEnabled = 1;
        
        for i=1:length(discreteEigenvalues)
            for j=1:length(discreteSpectrum)
                % Let's get the normalization parameters that we need to build the reference signal
                inft = INFT_v2(param.INFT);
                inft.normalizationParameters(const.c/Fc);
                
                % Autodetect A and B from the spectrum, theoretically we should be able to detect the phase as well
                %1-soliton
                A = imag(discreteEigenvalues(i)) * 2;
                B = A;
                
                % Generate an 1-soliton manually
                phase = angle(discreteSpectrum(i,j)) + pi/2;
                t0 = 1/A*log(abs(discreteSpectrum(i,j))/A);
                t = (-12:1/Fs/inft.Tn:12)*inft.Tn;
                x = A*sqrt(inft.Pn)*sech(B*(t/inft.Tn - t0))*exp(-1i*phase);
                sig = signal_interface(x, struct('Rs', inf, 'Fs', Fs, 'Fc', Fc));
                
                % Compute waveform with Darboux transform
                inft = DiscreteINFT_v1(param.INFT);
                sigDarb = inft.traverse(discreteEigenvalues(i), discreteSpectrum(i,j), Rs);
                
                % Compute NFT until all eignevalues are found. A little random.
                toten = 0;
                while toten < 0.9
                    nft_out = NFT_v8(param.NFT);
                    nft_out.traverse(sigDarb);
                    E = nft_out.results.E;
                    toten = (E.Ec + E.Ed)/E.ETot;
                end
                
                robolog('Discrete eigenvalues OUTPUT signal:');
                egDb = nft_out.discreteEigenvalues(); egDb
                robolog('Discrete spectrum');
                dsDb = nft_out.discreteSpectrum(); dsDb
                
                % Plot the results
                plotNFTConstellation_v1(egDb, dsDb, ...
                    'refEigenvalues', discreteEigenvalues(i), 'refSpectrum', discreteSpectrum(i,j));
                figure;
                inft.timeDomainSignal();
                subplot(1,2,1); hold on;
                plot(t, abs(x), 'r--');
                subplot(1,2,2); hold on;
                plot(t, angle(x)/pi, 'r--');
            end
        end
    case 6
        %% Test 6: Generate a 2-Solitions for a full 2x QPSK constellation
        
        nPoints = 2048; % Reduce this and you'll soon lose precision!
        Rs = Fs/nPoints;
        
        %INFT parameters
        param.INFT.Tn        = 0.03/Rs;                              %[km]
        param.INFT.gamma     = 1.27;                             %[/W/km]
        param.INFT.D         = 17;                            %[ps/(nm km)]
        param.INFT.method    = 'darboux';
        param.INFT.Fc        = Fc;
        param.INFT.nPoints   = nPoints;
        
        % Set the spectrum here
        discreteEigenvalues = 1i*[0.3 0.6];
        
        %NFT parameters
        param.NFT.Tn        = param.INFT.Tn;                              %[km]
        param.NFT.gamma     = param.INFT.gamma;                             %[/W/km]
        param.NFT.D         = param.INFT.D;                               %[ps/(nm km)]
        param.NFT.nPoints   = 2^10;
        param.NFT.method    = 'ABL';                          % FCD ABL N-ABL LP
        param.NFT.methodDiscreteSpectrum = 'ABL';
        param.NFT.computeDiscreteSpectrumEnabled = 1;
        param.NFT.complexPlaneSearchArea = 1.1 * (max(0.1, max(real(discreteEigenvalues))) + ...
            1i*max(imag(discreteEigenvalues)));
        param.NFT.mexEnabled = 1;
        
        %Right way
        %discreteSpectrum = -2*(discreteEigenvalues)*exp(1i*pi*[0:0.5:1.5]);
        %Buelow way
        discreteSpectrum = [exp(1i*pi*[0:0.5:1.5]); exp(1i*pi*[0:0.5:1.5])];
        
        discreteSpectrum(1,:) = discreteSpectrum(1,:)*exp(1i*0.25*pi);
        
        [a,b] = ndgrid(discreteSpectrum(1,:),discreteSpectrum(2,:));
        discreteSpectrum = [a(:), b(:)];
        %plotNFTConstellation_v1(discreteEigenvalues, discreteSpectrum);
        
        for j=1:size(discreteSpectrum,1)
            
            inft = DiscreteINFT_v1(param.INFT);
            sigDarb = inft.traverse(discreteEigenvalues, discreteSpectrum(j,:), Rs);
            
            % Compute NFT until all eignevalues are found. A little random.
            toten = 0;
            while toten < 0.9
                nft_out = NFT_v8(param.NFT);
                nft_out.traverse(sigDarb);
                E = nft_out.results.E;
                toten = (E.Ec + E.Ed)/E.ETot;
            end
            
            robolog('Discrete eigenvalues OUTPUT signal:');
            egDb = nft_out.discreteEigenvalues(); egDb
            robolog('Discrete spectrum');
            dsDb = nft_out.discreteSpectrum(); dsDb
            
            % Plot the results
            plotNFTConstellation_v1(egDb, dsDb, ...
                'refEigenvalues', discreteEigenvalues, 'refSpectrum', discreteSpectrum(j,:));
            figure;
            inft.timeDomainSignal();
        end
end
