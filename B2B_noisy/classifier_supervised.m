% This function classify and order eigenvalues found by NFT w.r.t. the
% original discrete spectrum

% WORKS IF THE ROUTINE FINDS MORE EIGS THAN THERE ARE

% @author Nicola De Renzis
% @date 16/02/2018

function [eigClass,ampClass] = classifier_supervised(eigNFT,ampNFT,eig,amp)

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

if N_est>1
    distance_eig = zeros(N_est,N);
    for i=1:N_est
        distance_eig(i,:) = abs(eig-eigNFT(i));
    end
    [~,min_locEigs] = min(distance_eig);
    spurious_locEigs = ~ismember((1:N_est),min_locEigs);
    eigClass = [eigNFT(min_locEigs),eigNFT(spurious_locEigs)];
    
    ampClass = [ampNFT(min_locEigs),ampNFT(spurious_locEigs)];
else
    eigClass = eigNFT;
    ampClass = ampNFT;
end

end

