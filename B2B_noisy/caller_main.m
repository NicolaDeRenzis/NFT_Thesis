clearall
close all

param_vec = [0.5*1i 1i 2*1i 4*1i];

for i=1:numel(param_vec)
    disp(numel(param_vec)-i)
    caller(param_vec(i))
end