function [q] = INFT_FastDarboux(eig,amp,h,N)
% author: Nicola De Renzis
% date: 19/01/2018
%
%
% This function implements the Fast INFT deployed at Delft Center
% [Sanders Wahls, Shrinivas Chimmalgi, Discrete Darboux based Fast Inverse
% Nonlinear Fourier Transform Algorithm for Multi-solitons]
%
% The algorithm is focused on adding discrete eigenvalues and
% generate the consequently waveform, which in this case (having only a discrete
% spectrum) is a multi-soliton wave.
%
% The algorithm is O(n*k) where n is the number of samples and k the number
% of eigenvalue. The precision is O(h^3).
%
%%% INPUT:
% eig:      Eigenvalues of the discrete spectrum. It can be a scalar or a
%           vector
% amp:      Spectral amplitude of eig in the discrete spectrum. It can be a
%           scalar or a vector. Must be consistent with eig.
% h:        time step calculated as follow: h = (T(N)-T(1))/N, where T is
%           time duration of the output signal and 2N+1 its number of
%           samples [s]
% N:        2N+1 is the number of samples of the signal in the time domain.
%
%%% OUTPUT
% q:        Output waveform in the time domain with time step 'h' and
%           length 2N+1 samples.

% check inputs
if size(amp)~=size(eig)
    error('The vectors of the eigenvalues and the amplitudes must have the same dimensions')
end
if isrow(eig)
    eig = eig.';
end
if isrow(amp)
    amp = amp.';
end

% iterate over all the eignvalue that you want to add.
K = length(eig);
[~,index] = sort(imag(eig),'descend');
eig = eig(index);
amp = amp(index);
z = exp(1i.*eig.*h);
b = amp;
% at n=0, where n indices the time of the waveform; -N<=n<=N
V0 = [1,0;0,1];
for j=1:K
    zj = z(j);  % current eigenvalue
    %bj = amp(j); % current amplitude? should I read from b or the original amplitudes?
    bj = b(j);
    za = zj;
    zb = 1/conj(zj);
    betaa = -bj;
    betab = 1/conj(bj);
    an(j) = (betab*zb-betaa*za)/(za*zb*(betaa*zb-betab*za));
    bn(j) = (za^2-zb^2)/(za*zb*(betaa*zb-betab*za));
    cn(j) = -conj(bn(j));
    dn(j) = conj(an(j));
    for k = j+1:K
        zk = z(k);
        b(k) = -(cn(j)-b(k)*(zk+dn(j)/zk))/(an(j)*zk+1/zk-b(k)*bn(j));
    end
    if j==1
        M0 = [zj^(-1)+an(j)*zj,bn(j);cn(j),zj+dn(j)^(-1)*zj^(-1)];
    else
        M0 = [zj^(-1)+an(1)*zj,bn(1);cn(1),zj+dn(1)^(-1)*zj^(-1)]*M0;
    end
end
V0p = M0*V0;
for n = 0:N
   Qp(N+1+n) = -bn(K)/dn(K);
   Rp = - conj(Qp(N+n));
   Lnp = [z(K)^(-1),Qp(N+n);Rp,z(K)];
   V0p = Lnp*V0p;
end
V0p = M0*V0;
for n=0:-N
   Qp(N+1+n) = bn(K);
   Rp = - conj(Qp(N+n));
   Lnp = 1/(1-Qp(N+n)*Rp).*[z(K)^(-1),-Qp(N+n);-Rp,z(K)];
   V0p = Lnp*V0p;
end

q = Qp./h;

end

