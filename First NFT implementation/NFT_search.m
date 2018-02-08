function [qc,lambda,qd,lambdaj] = NFT_search(q,t,varargin)
% author: Nicola De Renzis
% date: 15/01/2018
% update: 19/01/2018

% This script calculate the NFT of an input signal q(t) in time domain. The
% script calculate both the continuous and discrete domain, respectively
% qc(lambda) qd(lambdaj).
% Here are implemented several methods both for discrete and continuous
% spectra. The discrete is computed implementing a search method.
% Ref: [Information Trasmission using NFT Part II: numerical methods, Mansoor, Yousefi]

%%% INPUT
% q:    vector containing the time envelope of which the NFT will be
%       calculated
% t:    time axis
% varargin{1}: can assum different strings value for different methods:
%           'FW': Forward method
%           'CN': Crank-Nicolson method
%           'AL': Ablowitz-Ladik method
%           'LP': Layer-Peeling method (probably not suitable for discrete spectrum)
%%% OUTPUT
% qc:   continuous spectrum in NFT domain, vector
% lambda: eigenvalues of continuous spectrum, vector (Real mumber)
% qd:   discrete spectrum in NFT domain, vector
% lambdaj: discrete eigenvalues for discrete spectrum (Complex+ number)

if ~isempty(varargin)
    method = upper(varargin{1});
else
    method = 'FW'; % Forward
end

N = length(q);      % number of samples; index k
T1 = t(1);          % beginning of the support of the envelope q[k=0] [s]
T2 = t(end);        % end of the support of the envelope q[k=N] [s]
%T = t(end)-t(1);    % total duration of signal
eps = (T2-T1)/N;    % step in time [s]
realEigen_bound = 40; % edge of real eigenvalues grid
realEigen_step = 1e-1;% step of real eigenvalues grid

% define eigenvalues on real axis for qc
lambda = -realEigen_bound:realEigen_step:realEigen_bound;
%lambda = linspace(-8,8,8*1e2);
qc = zeros(1,length(lambda));
fprintf('Continuous Spectrum...\n')
% for l=1:length(lambda)
%     lcurr = lambda(l);  % current eigenvalue
%     v = computeEigenvectors(q,t,eps,lcurr,method);
%     qc(l) = v(2)/v(1)*exp(-2*1i*lcurr*t(N)); % store point for given eigenvalue (table Algorithm 1); eq (6a,6b,7)
% end

% compute conserved quantities eq. (20)

E1 = trapz(t,abs(q).^2);                                % energy of signal; better approx of integral.
E2 = 1/2/1i * trapz(t,q.*gradient(conj(q)));            % momentum of signal
E3 = -1/4 * trapz(t,abs(q).^4 - abs(gradient(q)).^2);   % Hamiltonian of signal

E1c = 1/pi * trapz(lambda,log(1+abs(qc).^2));           % continuous conserved quantity, k=1
E2c = 1/pi * trapz(lambda,lambda .* log(1+abs(qc).^2)); % continuous conserved quantity, k=2
E3c = 1/pi * trapz(lambda,lambda.^2 .* log(1+abs(qc).^2)); % continuous conserved quantity, k=3
if isnan(E1c)
    E1c = 0;
end

%initialize error
E = [E1-E1c,E2-E2c,E3-E3c]; % eq. (20)
err = norm(E(1));      % error of discrete eigenvalues; how much each eigenvalue improve the error
errSave(1) = err;
errPrev = Inf;      % save error of previous iteration
%
epsErr = 5*1e-1;        % error tolerance for all the eigenvalues (compared against err,errSave,errPrev)
epsDl = 1e-10;       % error toelarnce for each discrete eigenvalue to be considered as such
tol = 1e10*2.2204e-16; % tolerance for each eigenvalue if it has already been found
iterMax = 5*1e2;      % maximum number of iterations for each eigenvalue
eigenMax = 5*100;     % maximum number of eigenvalue (remember that the alg could converge on the same eigen)
timeoutMax = 100;    % number of iterations over which the overall error does not change
iter = 1;           % iteration counter
upperBoundD = 6;   % upper bound of complex plane (lower is 0)
upperBoundD = 1.2*max(real(q));   % upper bound of complex plane (lower is 0)
sideBoundD = 2;    % bound on the sides of the complex plane (symmetric)
alfa_i = 1;         % step parameter for Newton-Raphson method; set to 1 (19)
%
qd = zeros(1,eigenMax);
lambdaj = zeros(1,eigenMax);
NO_EIGENVALUES = true;
saveThisEigenv = false;
timeout_iter = 0;
counter_path = 1;

%while err>epsErr && iter<eigenMax && err<=errPrev && timeout_iter<=timeoutMax% each iteration will find a discrete eigenvalue.
while err>epsErr && iter<eigenMax && timeout_iter<=timeoutMax% each iteration will find a discrete eigenvalue.
    fprintf('Discrete Spectrum, iteration %d of %d\n',iter,eigenMax)
    i=1;            % inner while counter; mind here it starts at 1 instead of 0 as in the paper
    
    lambdaj_curr(i) = randFromRegion(-1/2,upperBoundD,-sideBoundD/2,sideBoundD/2); % initilize random starting eigenvalue
    path_stack(counter_path) = lambdaj_curr(i);counter_path=counter_path+1;
    condition = true;
    innerIter = 1;
    while condition % Newton-Raphson
        [v,dv] = computeEigenvectors(q,t,eps,lambdaj_curr(i),method);
        Dl = alfa_i * v(1)/(dv(1)+1i*t(N)*v(1));     % eq (19)
        lambdaj_curr(i+1) = lambdaj_curr(i)-Dl;
        path_stack(counter_path) = lambdaj_curr(i+1);counter_path=counter_path+1;
        figure(100);p(iter,innerIter)=plot(lambdaj_curr(1:i+1),'b-o','markersize',1);hold on;drawnow;title(sprintf('innerIter = %d, iter = %d',innerIter,iter));xlim([-sideBoundD/2-2 sideBoundD/2+2]);ylim([0-2 upperBoundD+2]);grid on;%pause
        if imag(lambdaj_curr(i+1))>upperBoundD || imag(lambdaj_curr(i+1))<=-2 || abs(real(lambdaj_curr(i+1)))>sideBoundD/2
            i=1;    % reset counter for new eigenvalue search
            lambdaj_curr(i) = randFromRegion(-1/2,upperBoundD,-sideBoundD/2,sideBoundD/2); % initilize random starting eigenvalue
            path_stack(counter_path) = lambdaj_curr(i);counter_path=counter_path+1;
        else
            i=i+1;
        end
        %
        %DlSave(innerIter)=abs(Dl);plot(1:innerIter,DlSave(1:innerIter));xlabel('iter');ylabel('error');drawnow
        if abs(Dl)<=epsDl  % if it enters here, then at least one time the algorithm has converged thanks to error decr.
            NO_EIGENVALUES = false;
            if imag(lambdaj_curr(i))>0
                set(p(iter,innerIter),'Color','red');figure(100);plot(lambdaj_curr(i),'go','markersize',3)
                saveThisEigenv = true; % keep track if the covergence is for an eigenvalue of our interest (IMAG > 0)
            end
        end
        innerIter = innerIter +1;
        condition = abs(Dl)>epsDl && innerIter<iterMax;
        TERM_DIST = abs(lambdaj_curr(i)-path_stack(1:counter_path-1-i));
        r = 1e-1; % radius for already-visited-point threshold
        %r=0;
        figure(100);h{iter,innerIter}=rectangle('Position',[[real(lambdaj_curr(i)),imag(lambdaj_curr(i))]-r 2*r 2*r],'Curvature',[1 1]);h{iter,innerIter-1}.EdgeColor = 'white';
        if any(TERM_DIST<r)  % check if the eig is close to already visited points
            condition = false;
            i=1;    % reset counter for new eigenvalue search
            lambdaj_curr(i) = randFromRegion(-1/2,upperBoundD,-sideBoundD/2,sideBoundD/2); % initilize random starting eigenvalue
            path_stack(counter_path) = lambdaj_curr(i);counter_path=counter_path+1;
            fprintf('      point already seen\n')
        end
    end % end of Netwon-Raphson
    
    if saveThisEigenv % now check that this eigenvalue is not already here
        ismember_eig = ismembertoli(lambdaj_curr(i),lambdaj,tol,'DataScale',1);
        if ~ismember_eig
            disp(lambdaj_curr(i))
            lambdaj(iter) = lambdaj_curr(i); % remember that lambdaj_curr will be overwritten maybe only partially
            qd(iter) = v(2)/(dv(1) + 1i*t(N)*v(1)) .* exp(-2*1i*lambdaj(iter)*t(N));
            bd(iter) = v(2).*exp(-1i*lambdaj(iter)*t(N)); % only the b(lambdaj)
            %update conditions
            Ed = [4*imag(lambdaj(iter)), 2*imag(lambdaj(iter)^2), 4/3*imag(lambdaj(iter)^3)];
            E = E - Ed;
            errPrev = err;
            err = norm(E(1));
            timeout_iter = 0;
        else
            timeout_iter = timeout_iter+1; % number of iterations in which the error does not change
        end
    else
        timeout_iter = timeout_iter+1; % number of iterations in which the error does not change
    end
    iter = iter+1;
    saveThisEigenv = false;
    
    errSave(iter)=err;figure(500);plot(1:iter,errSave(1:iter)); xlabel('iter');ylabel('error');drawnow
    
end

% take only eigenvalues further away than Matlab sensitivity
if NO_EIGENVALUES
    index = [];
    warning('No eigenvalues found');
else
    [~,index] = unique(lambdaj);
    index = index(2:end); % the first is always the 0, discard
    %[~,index] = find(lambdaj(index)>=0); % discard negative eigenvalues (in case)
end

% return final discrete spectrum
lambdaj = lambdaj(index);
qd = qd(index);
%%%%%%%%%%%%%%%%%%%%%
%qd=bd(index);
%%%%%%%%%%%%%%%%%%%%%


% subfunction to find the eigenvector through iterative procedure
    function [v,dv] = computeEigenvectors(func,t,step,currentEigenvalue,method)
        if strcmp(method,'LP') % Layer-Peeling
            % initial conditions, eq (10)
            a0 = 1;
            b0 = 0;
            da0 = 0;
            db0 = 0;
            D = @(n,L) sqrt(L^2 + abs(func(n))^2);
            x = @(n,L) (cos(D(n,L)*step)-1i*L/D(n,L)*sin(D(n,L)*step))*...
                exp(1i*L*(t(n)-t(n-1)));
            y = @(n,L) -conj(func(n))/D(n,L) .* sin(step*D(n,L))*exp(-1i*L*(t(n)+t(n-1)));
            dx = @(n,L) 1i*step*(1-L^2/(D(n,L))^2)*...
                (cos(D(n,L)*step)-sin(D(n,L)*step)/D(n,L)/step)*exp(1i*L*step);
            dy = @(n,L) -conj(func(n))*( L*step/(D(n,L))^2 *cos(D(n,L)*step) -...
                (L/(D(n,L))^3 + 1i * (t(n)+t(n-1))/D(n,L) )*sin(D(n,L)*step) )*...
                exp(-1i*L*(t(n)+t(n-1)));
            for k=2:N
                xk = x(k,currentEigenvalue);
                xk_c = conj(x(k,conj(currentEigenvalue)));
                yk = y(k,currentEigenvalue);
                yk_c = conj(y(k,conj(currentEigenvalue)));
                
                dxk = dx(k,currentEigenvalue);
                dxk_c = conj(dx(k,conj(currentEigenvalue)));
                dyk = dy(k,currentEigenvalue);
                dyk_c = conj(dy(k,conj(currentEigenvalue)));
                
                a = a0*xk - b0*yk_c;
                b = a0*yk + b0*xk_c;
                da = da0*xk + a0*dxk - (db0*yk_c+b0*dyk_c);
                db = da0*yk + a0*dyk + (db0*xk_c+b0*dxk_c);
                % update
                a0 = a;
                b0 = b;
                da0 = da;
                db0 = db;
            end
            
            v = [a*exp(-1i*currentEigenvalue*t(N)) , b*exp(1i*currentEigenvalue*t(N))];
            dv = [ da*exp(-1i*currentEigenvalue*t(N)) -1i*t(N)*v(1) ,...
                db*exp(1i*currentEigenvalue*t(N)) +1i*t(N)*v(2)];
        else
            
            %initial condition of eigenvector; eq (18a 18b)
            v0 = [1;0].*exp(-1i*currentEigenvalue*t(1));
            dv0 = [-1i*t(1);0].*exp(-1i*currentEigenvalue*t(1));
            for k=1:N-1
                if strcmp(method,'FW')
                    A = eye(2) + step.*[-1i*currentEigenvalue, func(k);-conj(func(k)), 1i*currentEigenvalue]; % eq (8)
                    dA = step.*[-1i,0;0,1i]; % eq. after (18b), derivative of A
                elseif strcmp(method,'CN') % Crank-Nicolson
                    A = (eye(2) - step/2.*[-1i*currentEigenvalue, func(k+1);-conj(func(k+1)) , 1i*currentEigenvalue])^(-1) *...
                        (eye(2) + step/2.*[-1i*currentEigenvalue, func(k);-conj(func(k)), 1i*currentEigenvalue]);
                    dA = step/2 .* (eye(2) - step/2.*[-1i*currentEigenvalue, func(k+1);-conj(func(k+1)), 1i*currentEigenvalue])^(-1)*...
                        diag([-1i,1i])*(eye(2) + A);
                else % Ablowitz-Ladik ('AL')
                    z = exp(-1i*currentEigenvalue*step);
                    A = [z,func(k)*step;-conj(func(k)*step),z^(-1)];
                    dA = step.*diag([-1i*z,1i*z^(-1)]);
                end
                % iterate to compute new eigenvector
                v = A*v0;
                dv = dA*v0 + A*dv0;
                v0 = v;
                dv0 = dv;
            end
        end
        
    end %subfunction

    function r = randFromRegion(Down,Up,Left,Right)
        % Complex number drawing from complex plane region
        %r = rand*(Right-Left)+Left + 1i*(rand*(Up-Down)+Down);
        r = rand*(Right-Left)+Left + 1i*(randn*sqrt(Up)*1.1).^2;
    end %subfunction
counter_path
end %function

