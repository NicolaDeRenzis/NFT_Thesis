% This function classify and order discrete amplitudes found by NFT w.r.t.
% needs to classify correctly the retreived eigenvalues


% @author Nicola De Renzis
% @date 26/02/2018

function [eigClass,ampClass,labels,flag_error] = createLabels(eigNFT,ampNFT,eig,amp)

% INPUT:
%   eigNFT      eigenvalues found by NFT
%   ampNFT      amplitudes found by NFT
%   eig         original eiganvalues
%   amp         original amplitudes

% OUTPUT:
%   eigClass    eigenvalues ordered
%   ampClass    amplitudes ordered

N = numel(eig); % number of eigs
N_est = numel(eigNFT); % number of estimated eigs
flag_error = 0;

if N_est>1
    imag_ampNFT = imag(ampNFT);
    for i=1:N % look for the closest ones to the real amplitudes
        diff = abs(imag_ampNFT - imag(amp(i)));
        [~,minLoc] = min(diff);
        imag_ampNFT_clean(i) = imag_ampNFT(minLoc);
        ampClass(i) = ampNFT(minLoc);
        eigClass(i) = eigNFT(minLoc);
    end
    [~,ind] = sort(imag_ampNFT_clean);
    ampClass = ampClass(ind);
    eigClass = eigClass(ind);
    labels = [1 2];
else
    eigClass = eigNFT;
    ampClass = ampNFT;
    labels = 1;
end

if sum(imag(eigClass)>1)~=1
    flag_error = 1;
end
end

