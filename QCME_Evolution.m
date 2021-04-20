%对于条件主方程，下一个时刻的密度矩阵和上一时刻相邻及自己的密度矩阵相关建立这样一个关系
clc;clear;close;format long
addpath 'G:\matlab\bin\QETLAB-0.9\QETLAB-0.9'
addpath G:\matlab\bin\QETLAB-0.9\QETLAB-0.9\helpers
rng('default')
para.n = 2;%密度矩阵的维度
para.s = 6;%一共有六个输出
rho_total = RandomDensityMatrix(para.n);
%先做一个简单的例子，假设密度矩阵总共被分为三个
a = rand(1);b = rand(1);c = 1-a-b;
rho_initial_zero = a.*rho_total;%n=0
rho_initial_one = b.*rho_total;%n=1
rho_initial_two = c.*rho_total;%n=2
rho.zero = rho_initial_zero;rho.one = rho_initial_one;rho.two = rho_initial_two;%分别对应第一，二，三个密度矩阵
disp('检验初始密度矩阵的有效性')

