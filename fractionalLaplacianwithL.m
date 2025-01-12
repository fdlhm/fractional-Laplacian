% fractionalLaplacianwithL.m
% If you use it, please cite the corresponding paper:
% Jorge Cayama, Francisco de la Hoz, Carlota Mar\'{\i}a Cuesta, Carlos Javier Garc\'{\i}a-Cervera,
% A fast convolution method for the fractional Laplacian in $\mathbb{R}$, (2025).
%
% This program approximates numerically the fractional Laplacian of
% $u(x)=\erf(x)$ ($\erf$ is defined in LaTeX by $\DeclareMathOperator{\erf}{erf}$).
% The variables N, r, a and L denote respectively $N$, $r$, $\alpha$ and $L$
clear
N=1048576;
r=8;
a=0.9;
L=2.1;
% sj is the vector containing the nodes $s_j$
sj=pi*((0:N-1)'+.5)/N;
% xj is the vector containing the nodes $x_j=\cot(s_j)$
xj=L*cot(sj);
tic
% uxj is the vector containing the values $u(x_j)=\erf(x_j)
uxj=erf(xj);
% Compute the FFT of an even extension of $\erf$
% u_ is the vector containing the values $\hat{u}(k)$
u_=fft([uxj;uxj(end:-1:1)]).*exp(-1i*pi*[0:N-1 -N:-1]'/(2*N));
u_(abs(u_)/(2*N)<eps)=0; % Apply the Krasny filter
% utilde_ is the vector containing the values $\hat{\tilde u}(k)$
utilde_=2*r*[u_(1:N);zeros((4*r-2)*N,1);u_(N+1:2*N)]...
    .*exp(1i*pi*[0:2*r*N-1 -2*r*N:-1]'/(4*r*N));
% ik is the vector containing $ik$, with
% $k\in\{0,\ldots,2rN-1\}\cup\{-2rN,\ldots,-1\}$
ik=1i*[0:2*r*N-1 -2*r*N:-1]';
% tildesn12 is the vector containing the values $\tilde s_{n+1/2}$
tildesn12=pi*((0:2*r*N-1)'+.5)/(2*r*N);
% ustildesn12 is the vector containing the numerical approximation of
% the values $u_s(\tilde s_{n+1/2})$
ustildesn12=ifft(ik.*utilde_);
% usstildesn12 is the vector containing the numerical approximation of
% the values $u_{ss}(\tilde s_{n+1/2})$
usstildesn12=ifft(ik.^2.*utilde_);
% F is the vector containing the numerical approximation of the values
% $f(\tilde s_{n+1/2}$
F=sin(tildesn12).*usstildesn12(1:2*r*N)...
    +2*cos(tildesn12).*ustildesn12(1:2*r*N);
% Invoke the Matlab function singularintegral, to obtain the vector I
% containing the numerical approximation of the values $I(s_j)$
I=singularintegral(N,r,a,1-a,F);
% fraclapnum and fraclap store respectively the numerical approximation
% and the exact value of the fractional Laplacian of $u(x)$
fraclapnum=(1/(L^a*2*gamma(2-a)*cos(pi*a/2)))*sin(sj).^(a-1).*I;
toc,tic % Elapsed time
fraclap=(2^(1+a)*gamma((1+a)/2)/pi)*xj.*hypergeom((1+a)/2,3/2,-xj.^2);
toc % Elapsed time
norm(fraclapnum-fraclap) % Error in discrete $L^2$ norm
norm(fraclapnum-fraclap,"inf") % Error in discrete $L^\infty$ norm
