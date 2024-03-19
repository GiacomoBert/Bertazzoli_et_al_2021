function R = corrm(A, B)
% Calculate the correlation matrix between two matrices
if size(A, 1) ~= size(B, 1)
    error('MATRIX DIMENSIONS MUST BE EQUAL!');
end
N1 = size(A, 2); N2 = size(B, 2);
Am = mean(A, 2); Bm = mean(B, 2);
Ar = A - repmat(Am, 1, N1); Br = B - repmat(Bm, 1, N2);
R = (diag(diag(Ar*Ar')))^(-0.5)*(Ar*Br')*(diag(diag(Br*Br')))^(-0.5);