clear all %#ok<CLALL>
close all

T1 = 0;
T2 = 1;
T = T2-T1;
A = 2;
di = 1e-6;
imagAxis = [0:di:A+1].*1i;
equation = exp(2*1i*(T2-T1).*sqrt(imagAxis.^2+abs(A).^2)) - (imagAxis + sqrt( imagAxis.^2 + A^2 )) ./ (imagAxis - sqrt( imagAxis.^2 + A^2 ));
figure(522);plot(imag(imagAxis),real(equation),'b-',imag(imagAxis),imag(equation),'b--');hold on

epsilon = di*2;
%e = epsilon +1i*epsilon;
[~,ind] = find(real(equation)>-epsilon & real(equation)<epsilon & imag(equation)>-epsilon & imag(equation)<epsilon );
k=1; % compare position offset
j=1; % current position
i=1; % free position pointer;
while j+k<=length(ind)
    if ind(j)+k==ind(j+k)
        if j+k==length(ind)
            ind(i) = ind(j+ceil((k-1)/2));
        end
        k=k+1;
    else
        ind(i) = ind(j+ceil((k-1)/2));
        i=i+1;
        if j+k==length(ind)
            ind(i) = ind(end);
        end
        j=j+k;
        k=1;
    end
end
ind = ind(1:i);
figure(522);plot(imag(imagAxis(ind)),zeros(1,length(ind)),'ro');hold off
eig = imagAxis(ind(1:end-1));

