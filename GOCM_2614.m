%在复数空间产生4类Kraus算符，每类算符一共有六个
clc;clear;
clc;clear;close;
para.n = 2;para.s = 6; para.w = 1;para.class = 4;
m = para.n * para.s * para.w * para.class;
seed = 1;rng(seed);
coL1 = rand(m, 1);coL2 = rand(m, 1);
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
coL2 = coL2 - dot(coL1, coL2)/norm(coL1)^2.*coL1;

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
coL3 = coL3 - dot(coL1,coL3)/norm(coL1)^2.*coL1-dot(coL2,coL3)/norm(coL2)^2.*coL2;

coL4 = rand(m, 1);
RREF = rref([coL1';coL2';coL3';coL4']);
for i = 1:m
    if RREF(3,i) == 0
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
coL4 = coL4 - dot(coL1,coL4)/norm(coL1)^2.*coL1 - dot(coL2,coL4)/norm(coL2)^2.*coL2 - dot(coL3,coL4)/norm(coL3)^2.*coL3;

CoL1 = coL1 + 1i* coL2;
CoL2 = coL3 + 1i* coL4;

CoL1 = CoL1./norm(CoL1);
CoL2 = CoL2./norm(CoL2);
Kappa = [CoL1,CoL2];
Kappa_real = real(Kappa);
Kappa_imag = imag(Kappa);
save Kappa_real.txt -ascii Kappa_real
save Kappa_imag.txt -ascii Kappa_imag