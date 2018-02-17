% variance checker

clear all %#ok<CLALL>
close all

% remember, OSNR decreasing

tmp{1} = load('variance_0_5i.mat');
tmp{2} = load('variance_1_0i.mat');
tmp{3} = load('variance_2_0i.mat');
tmp{4} = load('variance_4_0i.mat');

for i=2:4
    var_ratio_eig(i-1,:) = tmp{i}.var_errorEigs(:,1) ./ tmp{i-1}.var_errorEigs(:,1);
    var_ratio_amp(i-1,:) = tmp{i}.var_errorAmpl(:,1) ./ tmp{i-1}.var_errorAmpl(:,1);
end