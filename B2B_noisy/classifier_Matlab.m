clear all %#ok<CLALL>
close all

load store_1000realiz_1.5i_2.mat

osnrLevel = 8;
eigsTot = store{osnrLevel}.eigs;
amplTot = store{osnrLevel}.ampl;

eigReal = [real(eigsTot(:,1));real(eigsTot(:,2))];
eigImag = [imag(eigsTot(:,1));imag(eigsTot(:,2))];
ampReal = [real(amplTot(:,1));real(amplTot(:,2))];
ampImag = [imag(amplTot(:,1));imag(amplTot(:,2))];
label = [zeros(length(eigsTot),1);ones(length(eigsTot),1)];

data = [eigReal,eigImag,ampReal,ampImag,label];