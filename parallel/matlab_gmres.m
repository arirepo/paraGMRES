%script for doing gmres 
clear all
clc

epsil = 1.e-10;

%reading the command line arguments
filename = '../matrices/dense6x6.mtx';
	%openning the input file
fid = fopen(filename, 'r');
	%reading the header
tmp = fscanf(fid, '%d %d %d\n', 3);
nrows = tmp(1);
nnz = tmp(3);

	%generates sparse matrix
A = sparse(nrows, nrows);
	%scanning through the file and constructing sparse matrix on fly     
for i = 1:nnz 
  tmp = fscanf(fid, '%d %d %lg\n', 3);
  A(tmp(1), tmp(2)) = tmp(3);
end
%close the input file
fclose(fid);

	%generates rhs vector for debugging
    %reading the command line arguments
filename = '../matrices/A_rhs.mtx';
	%openning the input file
fid = fopen(filename, 'r');
	%reading the header
tmp = fscanf(fid, '%d %d\n', 2);

	%generates rhs
b = sparse(nrows, 1);
	%scanning through the file and constructing sparse matrix on fly     
for i = 1:nrows 
  tmp = fscanf(fid, '%lg\n', 1);
  b(i, 1) = tmp(1);
end
%close the input file
fclose(fid);

b = sparse(nrows,1);
for i = 1:nrows
   b(i,1) = 1.;
end

%number of gmres iterations
n = nrows;
num = nrows;
v=zeros(nrows,num);
x0 = ones(nrows,1);
%x0 = x0/norm(x0);
r0 = b- A*x0;

v(:,1) = r0/norm(r0);
g(1) = norm(r0);

for k = 1:num
    g(k+1) = 0.;

    uk = A*v(:,k);
    for j = 1:k
        h(j,k) = (v(:,j)') * uk;
        uk = uk - h(j,k) * v(:,j);
    end
    
    h(k+1,k) = norm(uk);
    
    v(:,k+1) = uk/h(k+1,k);
   
    for j = 1:(k-1)
        delta = h(j,k);
        h(j,k) = c(j)*delta + s(j)*h(j+1,k);
        h(j+1,k) = -s(j)*delta + c(j) * h(j+1,k);
    end
    gamma = sqrt(h(k,k)^2 + h(k+1,k)^2);
    c(k) = h(k,k) / gamma;
    s(k) = h(k+1,k) / gamma;
    h(k,k) = gamma;
    h(k+1,k) = 0.;
    delta = g(k);
    g(k) = c(k) * delta + s(k) * g(k+1);
    g(k+1) = -s(k) * delta + c(k) * g(k+1)
    resi(k) = abs(g(k+1)); 
    if( resi(k) < epsil) 
       alpha = h(1:end-1,1:end)\(g(1:(end-1))');
       zk = zeros(nrows,1);
       for j = 1:k
           zk = zk + alpha(j)*v(:,j);
       end
       x(:,k+1) = x0 + zk;
       x(:,k+1)
       hold on;
       semilogy(resi);
       xlabel('Iterations');
       ylabel('Residuals');
       grid on;
       return 
    end
 
end


%format short e
%full(h)
%axis equal
%pause
    