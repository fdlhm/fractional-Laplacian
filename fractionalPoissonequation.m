% fractionalpoissonequation.m
% fractionalLaplacianfunction.m
% If you use it, please cite the corresponding paper:
% Jorge Cayama, Francisco de la Hoz, Carlota Mar\'{\i}a Cuesta, Carlos Javier Garc\'{\i}a-Cervera,
% A fast convolution method for the fractional Laplacian in $\mathbb{R}$, (2025).
%
% This program solves the fractional Poisson equation, i.e., given $f(x)$,
% it finds $u(x)$, such that $(-\Delta^{\alpha/2})u(x)=f(x)$, and $u(-\infty)=u_0$.
% In this example, $f(x)=\erf(x)$ ($\erf$ is defined in LaTeX by $\DeclareMathOperator{\erf}{erf}$).
% The variables N, r, a and L denote respectively $N$, $r$, $\alpha$ and $L$
clear
N=128;
r=256;
a=0.1;
L=2.1;
% sj is the vector containing the nodes $s_j$
sj=pi*((0:N-1)'+.5)/N;
% xj is the vector containing the nodes $x_j=\cot(s_j)$
xj=L*cot(sj);
tic
% uxj is the vector containing the exact values $u(x_j)=\erf(x_j)$
uxj=erf(xj);
u0=-1;
% fxj is the vector containing the values $f(x_j)$
fxj=(2^(1+a)*gamma((1+a)/2)/pi)*xj.*hypergeom((1+a)/2,3/2,-xj.^2);
% Ma is the operational matrix $M_\alpha$
Ma=zeros(N,N);
% I is the identity matrix of order N, whose $k$th column is precisely $\mathbf e_k$
I=eye(N);
% b is the row vector containing the values of $u_k(-\infty)$
b=zeros(1,N);
for k=1:N
    uk=I(:,k); % $u_k(x)$
    [fraclapuk,ukminusinfinity]=fractionalLaplacianfunction(I(:,k),r,a,L); % Create the kth column of Ma;
    Ma(:,k)=fraclapuk;
    b(k)=ukminusinfinity;
end
% uxj is the vector containing the numerical approximation of the values $u(x_j)$
uxjnum=[Ma;b]\[fxj;u0];
norm(uxjnum-uxj,inf) % Error in discrete $L^2$ norm
norm(uxjnum-uxj,2) % Error in discrete $L^\infty$ norm
toc
