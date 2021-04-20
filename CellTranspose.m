function y = CellTranspose(A)
len = length(A);
y = cell(len,1);
for i = 1:1:len
    y{i,1} = A{i}';
end
end