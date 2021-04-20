%验证密度矩阵循环迭代产生结果大小
clc;clear ;close
addpath 'G:\matlab\bin\QETLAB-0.9\QETLAB-0.9'
addpath G:\matlab\bin\QETLAB-0.9\QETLAB-0.9\helpers
rho0 = RandomDensityMatrix(2,'brues');
Kappa_real = load("Kappa_real261.txt");
Kappa_imag = load("Kappa_imag261.txt");
Kappa = Kappa_real + 1i * Kappa_imag;
seq = load("sequence.txt");
train_seq = seq(1,:);%取第一行带入计算
prob = [];
%%
%使用迭代法来求似然函数
for i = 1:1:length(train_seq)
    j = train_seq(i);
    K = Kappa(2*j-1:2*j,:);
    rho1 = K*rho0*K';
    P = trace(rho1);
    prob(i) = P;
    rho1 = rho1./P;
    rho0 = rho1;  
end
plot(1:1:length(train_seq),prob)

%%
%用公式法来求解
for i = 1:1:length(train_seq)
    j = train_seq(i);
    K = Kappa(2*j-1:2*j,:); 
    rho1 = K*rho0*K';
    rho0 = rho1;
end
disp(trace(rho0))