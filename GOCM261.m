%在复数空间中产生2*2的Kraus算符，一共有6个
clc;clear;close;format long
seed = 1;
rng(seed)
Para.n = 4;%实数空间的维数
Para.s = 6;%Kraus算符的个数
Para.w = 1;%每个输出对应的Kraus算符的个数
m = Para.n/2 * Para.s * Para.w;
coL1 = rand(m,1);coL1 = coL1./norm(coL1);%第一列
coL2 = rand(m,1);%随机产生第一列，判断与第一列线性无关，然后在施密特正交化
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
coL3 = rand(m,1);%第三列
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
coL3 = coL3 - dot(coL1,coL3)/norm(coL1)^2.*coL1-dot(coL2,coL3)./norm(coL2)^2.*coL2;coL3 = coL3./norm(coL3);
coL4 = rand(m,1);%第四列
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
coL4 = coL4 - dot(coL1,coL4)/norm(coL1)^2.*coL1-dot(coL2,coL4)/norm(coL2)^2.*coL2 - dot(coL3,coL4)/norm(coL3)^2.*coL3;
coL4 = coL4./norm(coL4);
%组合实数空间中的向量，放在复数空间中
CoL1 = coL1 + 1i*coL2;CoL1 = CoL1./norm(CoL1);
CoL2 = coL3 + 1i*coL4;CoL2 = CoL2./norm(CoL2);
Kappa = [CoL1,CoL2];
%验证Kraus算符的性质
Kraus_Operator = [];
Kraus = zeros(Para.n /2, Para.n / 2);
j = 1;
for k = 1:2:m
    Kraus = Kappa(k:k+1,:)'*Kappa(k:k+1,:) + Kraus;
    Kraus_Operator{j} = Kappa(k:k+1,:);
    j = j+1;
    disp(Kappa(k:k+1,:))
    
end
fprintf('The sum of all Karus dagger multi Kraus matrix(complex)\n')
disp(Kraus)
Kappa_real = real(Kappa);
Kappa_imag = imag(Kappa);
%save Kappa_real261.txt -ascii Kappa_real
%save Kappa_imag261.txt -ascii Kappa_imag



%%
seq = 1 + 5.*rand(1,10);
seq = round(seq);
addpath 'G:\matlab\bin\QETLAB-0.9\QETLAB-0.9'
addpath G:\matlab\bin\QETLAB-0.9\QETLAB-0.9\helpers
rho = RandomDensityMatrix(Para.n/2);
temp = eye(2);
%计算like-hood function
for i = 1:length(seq)
    temp = Kraus_Operator{seq(i)}*temp;
end
L = -log(real(trace(temp*rho*temp')));
fprintf('The value of likelihood function :\n')
disp(L)
save seq.txt -ascii seq