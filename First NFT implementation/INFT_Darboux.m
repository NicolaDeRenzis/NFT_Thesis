function q = INFT_Darboux(eig,amp,t)
% author: Nicola De Renzis
% date: 25/01/2018
%
%
% This function implements the INFT based on Darboux transform
% [Control and Detection of discrete spectra amplitudes in NL Fourier Spectrum]
% Algorithm 2
%
% The algorithm is focused on adding discrete eigenvalues and
% generate the consequently waveform, which in this case (having only a discrete
% spectrum) is a multi-soliton wave.
%
% The algorithm is O(N^2) where N is the number of eigenvalues.
%
%%% INPUT:
% eig:      Eigenvalues of the discrete spectrum. It can be a scalar or a
%           vector
% amp:      Discrete nonlinear spectrum b(lambdai)/da(lambdai). It can be a
%           scalar or a vector. Must be consistent with eig.
% t:        time axis of the output waveform.
%
%%% OUTPUT
% q:        Output waveform in the time domain with time step 'h' and
%           length of t samples.

% check inputs
if size(amp)~=size(eig)
    error('The vectors of the eigenvalues and the amplitudes must have the same dimensions')
end
if ~isrow(eig)
    eig = eig.';
end
if ~isrow(amp)
    amp = amp.';
end
if ~isrow(t)
    t = t.';
end

% iterate over all the eignvalue that you want to add.
N = length(eig);
n = length(t); % number of samples
q = zeros(1,n);
TERM1 = amp./(eig-conj(eig)); % put this outside the for loop because is invariant (need a matrix)
%%{
rho = zeros(N,n);
for i=1:N
    lk = [eig(1:i-1),eig(i+1:end)];
    TERM2 = (eig(i)-lk)./(eig(i)-conj(lk));
    rho(i,:) =  TERM1(i).*prod(TERM2).*exp(2*1i*eig(i).*t); % this is rho^(0)
    %rho(i,:) =  amp(i).*exp(2*1i*eig(i).*t); % WORKS IF YOU INPUT bk. this is rho^(0)
end
old_rho = rho; % old_rho --> rho^(i-1), rho --> rho^(i)
for i=1:N
    %rho = zeros(N,n); % shouldn't be necessary
    li = eig(i);
    curr_rho = old_rho(i,:);
    TERM3 = 1+abs(curr_rho).^2;
    q = q+2*1i*(li-conj(li)).*conj(curr_rho)./TERM3;
    for k=i+1:N
        rho(k,:) = ( (eig(k)-li).*old_rho(k,:)+(li-conj(li))./TERM3 .* (old_rho(k,:)-curr_rho) ) ./...
            ( eig(k)-conj(li) - (li-conj(li))./TERM3.* (1+conj(curr_rho).*old_rho(k,:) )  ) ;
    end
    old_rho = rho;
end
%}

%{
rho = zeros(N,n,N+1);
for i=1:N
    lk = [eig(1:i-1),eig(i+1:end)];
    TERM2 = (eig(i)-lk)./(eig(i)-conj(lk));
    rho(i,:,1) =  TERM1(i).*prod(TERM2).*exp(2*1i*eig(i).*t); % this is rho^(0)
end
%old_rho = rho; % old_rho --> rho^(i-1), rho --> rho^(i)
for i=1:N
    li = eig(i);
    curr_rho = rho(i,:,i);
    TERM3 = 1+abs(curr_rho).^2;
    q = q+2*1i*(li-conj(li)).*conj(curr_rho)./TERM3;
    for k=i+1:N
        rho(k,:,i+1) = ( (eig(k)-li).*rho(k,:,i)+(li-conj(li))./TERM3 .* (rho(k,:,i)-curr_rho) ) ./...
            ( eig(k)-conj(li) - (li-conj(li))./TERM3.* (1+conj(curr_rho).*rho(k,:,i) )  ) ;
    end
end
%}

%{
for i=1:N
    A = 1;
    %     product=1;
    %     for k=1:N
    %        if k~=i
    %            product = product*(eig(i)-eig(k))/(eig(i)-conj(eig(k)));
    %        end
    %     end
    %     B =  -TERM1(i).*product;
    lk = [eig(1:i-1),eig(i+1:end)];
    TERM2 = (eig(i)-lk)./(eig(i)-conj(lk));
    B =  -TERM1(i).*prod(TERM2);
    %B = -amp(i);
    v1(i,:) = A*exp(-1i*eig(i).*t);
    v2(i,:) = B*exp(+1i*eig(i).*t);
end
cond = false;
for i=1:N
    li = eig(i);
    abs_val = abs(v1(i,:)).^2 + abs(v2(i,:)).^2;
    q = q + 2*1i*(conj(li)-li).*(conj(v2(i,:) .* v1(i,:)))./abs_val;
    
    for k=i+1:N
        lk = eig(k);
        v1_new(k,:) = (lk-conj(li) - (li-conj(li)).*abs(v1(i,:)).^2./abs_val).*v1(k,:) -...
            (li-conj(li)).*conj(v2(i,:)).*v1(i,:)./abs_val.*v2(k,:);
        v2_new(k,:) = -(li-conj(li)).*v2(i,:).*conj(v1(i,:))./abs_val.*v1(k,:) +...
            (lk-li + (li-conj(li)).*abs(v1(i,:)).^2./abs_val).*v2(k,:);
        cond = true;
    end
    if cond
        v1 = v1_new;
        v2 = v2_new;
    end
end
%}

%% SIMONE
%{

L=n;
lambdas = eig;
discreteSpectrum = amp;
normByColumn = @(x) sqrt(sum(abs(x).^2,1));
obj.setNFTParameterB = 0;
qs = zeros(1, L); % Array of row vectors. The n-th vector is the time signal corresponding
% to the spectrum containing the first n eigenvalues
v = zeros(2, L, N); % Column vector with 2 elements x L (one per each element of qs) x N eigenvalues

% Init the eigenvectors for each eigenvalue [\ref Aref2016]
for i=1:N
    A = 1;
    if obj.setNFTParameterB
        B = -discreteSpectrum(i);
    else
        tempProduct = 1;
        for k=1:N
            if k ~= i
                tempProduct = tempProduct*(lambdas(i) - lambdas(k))/(lambdas(i) - conj(lambdas(k)));
            end
        end
        B = -discreteSpectrum(i)/(lambdas(i) - conj(lambdas(i)))*tempProduct;
    end
    v(1, :, i) = A*exp(-1i*lambdas(i)*t);
    v(2, :, i) = B*exp(1i*lambdas(i)*t);
end

% The block diagram at mansoor Part III p.11 explains more clearly what's happening here.

%For each eigenvalue..
for k=1:N
    qs = qs + 2*1i*(conj(lambdas(k)) - lambdas(k)) ...
        .*v(1,:,k).*conj(v(2,:,k))./normByColumn(v(:,:,k)).^2;
    
    vnew = zeros(2, L, N);
    %Update the eigenvalues.. (not updated for K = N)
    for j = k+1:N
        % The eigenvalues are update using the equation in the table in Mansoor Part III, p.12
        % This equations are just an algebraic semplification of the matrix equation:
        % u = (lambda*I - S Gamma S^-1) * v      [eq.12]
        vnew(1,:,j) = 1./normByColumn(v(:,:,k)).^2 .* (( ...
            (lambdas(j) - lambdas(k)) .* normByColumn(v(1,:,k)).^2 + ...
            (lambdas(j) - conj(lambdas(k))) .* normByColumn(v(2,:,k)).^2 ...
            ) .* v(1,:,j) + ...
            (conj(lambdas(k)) - lambdas(k)).*v(1,:,k).*conj(v(2,:,k)).*v(2,:,j));
        vnew(2,:,j) = 1./normByColumn(v(:,:,k)).^2 .* (( ...
            (lambdas(j) - conj(lambdas(k))) .* normByColumn(v(1,:,k)).^2 + ...
            (lambdas(j) - lambdas(k)) .* normByColumn(v(2,:,k)).^2 ...
            ) .* v(2,:,j) + ...
            (conj(lambdas(k)) - lambdas(k)).*conj(v(1,:,k)).*v(2,:,k).*v(1,:,j));
    end
    v = vnew;
end
q = qs.'; % Return a column vector
%}




end

