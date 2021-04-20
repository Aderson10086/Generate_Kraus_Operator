%产生一个2*2的四个复数空间的随机的Kraus算符
%具体的产生方法参看代码‘GenerateOrthonormalColumnsMaxtrix.m’文件
clc;clear;
n = 4;s=2;w=1;m = n*s*w;
rng('default')
coL1 = rand(m,1);coL1 = coL1./norm(coL1);
coL2 = rand(m,1);
RREF = rref([coL1';coL2']);
for i = 1:m
    if RREF(2,i) == 0
        if i == m
            disp('The first two columns are linear-dependent to each other')
            disp('please retry')
        end
        continue
    else
        disp('The first two columns are linear-independent to each other')
        break
    end
end
coL2 = coL2 - dot(coL1,coL2)/norm(coL1)^2.*coL1;coL2 = coL2./norm(coL2);
coL3 = rand(m,1);
RREF = rref([coL1';coL2';coL3']);
for i = 1:m
    if RREF(3,i) == 0
        if i == m
            disp('The first three columns are linear-dependent to each other')
            disp('please retry')
        end
        continue
    else
        disp('The first three columns are linear-independent to each other')
        break
    end
end
coL3 = coL3 - dot(coL3,coL2)/norm(coL2)^2.*coL2-dot(coL3,coL1)/norm(coL1)^2.*coL1;
coL3 = coL3/norm(coL3);
coL4 = rand(m,1);
RREF = rref([coL1';coL2';coL3';coL4']);
for i = 1:m
    if RREF(4,i) == 0
        if i == m
            disp('The first four columns are linear-dependent to each other')
            disp('please retry')
        end
        continue
    else
        disp('The first four columns are linear-independent to each other')
        break
    end
end
coL4 = coL4 - dot(coL4,coL3)/norm(coL3)^2.*coL3-dot(coL4,coL2)/norm(coL2)^2.*coL2-dot(coL4,coL1)/norm(coL1)^2.*coL1;
coL4 = coL4./norm(coL4);
j = sqrt(-1);
CoL1 = coL1 + j * coL2;CoL1 = CoL1./norm(CoL1);
CoL2 = coL3 + j * coL4;CoL2 = CoL2./norm(CoL2);
Kappa = [CoL1,CoL2];
fprintf('vertically statck Kraus operator\n')
disp(Kappa)
Kraus = zeros(n/2,n/2);
for k = 1:2:n*s*w
    Kraus = Kappa(k:k+1,:)'*Kappa(k:k+1,:) + Kraus;
end
fprintf('The sum of all Karus dagger multi Kraus matrix(complex)\n')
disp(Kraus)
Kappa_real = real(Kappa);
Kappa_imag = imag(Kappa);
% save Kappa_real.txt -ascii Kappa_real
% save Kappa_imag.txt -ascii Kappa_imag
addpath 'G:\matlab\bin\QETLAB-0.9\QETLAB-0.9'
addpath G:\matlab\bin\QETLAB-0.9\QETLAB-0.9\helpers
rho = RandomDensityMatrix(s);
rho_real = real(rho);
% save rho_real.txt -ascii rho_real
% rho_imag = imag(rho);
% save rho_imag.txt -ascii rho_imag