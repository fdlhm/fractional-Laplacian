% fractionalLaplacianexamplewithoutL.m
% If you use it, please cite the corresponding paper:
% Jorge Cayama, Francisco de la Hoz, Carlota Mar\'{\i}a Cuesta, Carlos Javier Garc\'{\i}a-Cervera,
% A fast convolution method for the fractional Laplacian in $\mathbb{R}$, (2025).
%
% This program approximates numerically the fractional Laplacian of 
% $u(x)=(ix-1)/(ix+1)$. We work with $u(s)\equiv u(\cot(s))=e^{2is}$
% The variables N, r and a denote respectively $N$, $r$ and $\alpha$
clear
tic
N=10000019;
r=1;
a=1.3;
% Define the function
% $f(s)=\sin(s)u_{ss}(s)+2\cos(s)u_{s}(s)=(-4\sin(s)+4i\cos(s))e^{2is}$
f=@(s)(-4*sin(s)+4i*cos(s)).*exp(2i*s);
% F is the vector containing the values $f(\tilde s_{n+1/2}$
F=f(pi*((0:2*r*N-1)'+1/2)/(2*r*N));
% Invoke the Matlab function singularintegral, to obtain the vector I
% containing the numerical approximation of the values $I(s_j)$
I=singularintegral(N,r,a,1-a,F); % Invoke the function singularintegral
% sj is the vector containing the nodes $s_j$
sj=pi*((0:N-1)'+.5)/N;
% fraclapnum and fraclap store respectively the numerical approximation
% and the exact value of the fractional Laplacian of $u(x)$
fraclapnum=(1/(2*gamma(2-a)*cos(pi*a/2)))*sin(sj).^(a-1).*I;
fraclap=-2*gamma(1+a)./(1i*cot(sj)+1).^(1+a);
norm(fraclapnum-fraclap) % Error in discrete $L^2$ norm
norm(fraclapnum-fraclap,"inf") % Error in discrete $L^\infty$ norm
toc % Elapsed time
