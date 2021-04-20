function y = CellMultiDagger(A,B)
y = zeros(2,2);C = CellTranspose(A);
for i = 1:1:length(A)
    y = y + A{i}*B{i}*C{i};
end