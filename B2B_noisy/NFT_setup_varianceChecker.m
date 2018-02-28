% variance checker

clear all %#ok<CLALL>
close all

% remember, OSNR decreasing

% tmp{1} = load('store_2000realiz_0.5i.mat');
% tmp{2} = load('store_2000realiz_1.0i.mat');
% tmp{3} = load('store_2000realiz_2.0i.mat');
% tmp{4} = load('store_2000realiz_4.0i.mat');
tmp{1} = load('store_100realiz_2.0i.mat');
%%
%{
for i=2:4
    var_ratio_eig(i-1,:) = tmp{i}.var_errorEigs(:,1) ./ tmp{i-1}.var_errorEigs(:,1);
    var_ratio_amp(i-1,:) = tmp{i}.var_errorAmpl(:,1) ./ tmp{i-1}.var_errorAmpl(:,1);
end
%}
%%
ampl_residues = [];
eigs_residues = [];
for i=1:1 % indices the eigs imaginary part amplitude
    local = tmp{i};
    for k=1:length(local.store) % indices the osnr level
        eigs = local.store{k}.eigs;
        eigs_orig = local.store{k}.eigs_original;
        firstSpurious = length(eigs_orig)+1;
        ampl = local.store{k}.ampl;
        ampl_orig = local.store{k}.ampl_original;
        index_spurious = find(eigs(:,firstSpurious)); % if any, find the non-zero elements
        if ~isempty(index_spurious)
            temp = ampl(:,1:firstSpurious)-[repmat(ampl_orig,length(ampl),1),zeros(length(ampl),1)];
            ampl_residues = [temp(index_spurious,:);ampl_residues];
            temp = eigs(:,1:firstSpurious)-[repmat(eigs_orig,length(eigs),1),zeros(length(eigs),1)];
            eigs_residues = [temp(index_spurious,:);eigs_residues];
        end
    end
end

residue_ampl_diff = abs(ampl_residues(:,2))-abs(ampl_residues(:,1));
negative_residue_ampl_diff = find(residue_ampl_diff<=0);

residue_eigs_diff = abs(eigs_residues(:,2))-abs(eigs_residues(:,1));
negative_residue_eigs_diff = find(residue_eigs_diff<=0);

% the second column is classified as spuripous eigen, so it should have a
% higher residue compared to the first one, which is classified as correct
% eigenvalue. If it is not the case then it's a problem
figure
plot(residue_ampl_diff)
if ~isempty(negative_residue_ampl_diff)
    hold on
    plot(negative_residue_ampl_diff,residue_ampl_diff(negative_residue_ampl_diff),'o')
end
plot(residue_eigs_diff)
if ~isempty(negative_residue_ampl_diff)
    %hold on
    plot(negative_residue_ampl_diff,residue_eigs_diff(negative_residue_ampl_diff),'o')
end
hold off
