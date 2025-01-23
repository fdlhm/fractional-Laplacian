% fractionalLaplacianfunction.m
% If you use it, please cite the corresponding paper:
% Jorge Cayama, Francisco de la Hoz, Carlota Mar\'{\i}a Cuesta, Carlos Javier Garc\'{\i}a-Cervera,
% A fast convolution method for the fractional Laplacian in $\mathbb{R}$, (2025).
%
% This program computes the fractional Laplacian of a function $u(x)$
% whose values $u(x_j)$, stored in the vector uxj, are known.
% The variables N, r, a and L denote respectively $N$, $r$, $\alpha$ and $L$
function fraclapnum=fractionalLaplacianfunction(uxj,r,a,L)
N=length(uxj);
% sj is the vector containing the nodes $s_j$
sj=pi*((0:N-1)'+.5)/N;
% u_ is the vector containing the values $\hat{u}(k)$
u_=fft([uxj;uxj(end:-1:1)]).*exp(-1i*pi*[0:N-1 -N:-1]'/(2*N));
u_(abs(u_)/(2*N)<eps)=0; % Apply the Krasny filter
% utilde_ is the vector containing the values $\hat{\tilde u}(k)$
utilde_=2*r*[u_(1:N);zeros((4*r-2)*N,1);u_(N+1:2*N)].*exp(1i*pi*[0:2*r*N-1 -2*r*N:-1]'/(4*r*N));
% ik is the vector containing $ik$, with $k\in\{0,\ldots,2rN-1\}\cup\{-2rN,\ldots,-1\}$
ik=1i*[0:2*r*N-1 -2*r*N:-1]';
tildesn12=pi*((0:2*r*N-1)'+.5)/(2*r*N); % $\tilde s_{n+1/2}$
Ustildesn12=ifft(ik.*utilde_); % $U_s(\tilde s_{n+1/2})$
Usstildesn12=ifft(ik.^2.*utilde_); % $U_{ss}(\tilde s_{n+1/2})$
% F is the vector storing the numerical approximations of
% $f(s)=\sin(s)u_{ss}(s)+2\cos(s)u_{s}(s)$ at $s_j$
F=sin(tildesn12).*Usstildesn12(1:2*r*N)+2*cos(tildesn12).*Ustildesn12(1:2*r*N);
I=singularintegral(N,r,a,1-a,F); % Invoke the function singularintegral
% fraclapnum stores the numerical approximation of the fractional Laplacian of $u(x)$
fraclapnum=(1/(L^a*2*gamma(2-a)*cos(pi*a/2)))*sin(sj).^(a-1).*I;
