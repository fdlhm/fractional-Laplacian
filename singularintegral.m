% singularintegral.m
% If you use it, please cite the corresponding paper:
% Jorge Cayama, Francisco de la Hoz, Carlota Mar\'{\i}a Cuesta, Carlos Javier Garc\'{\i}a-Cervera,
% A fast convolution method for the fractional Laplacian in $\mathbb{R}$, (2025).
%
% This function approximates numerically
% I(s_j)=\int_0^\pi\sin^\beta(\eta)|\sin(\eta-s_j)|^\gamma f(\eta)d\eta
% The parameters N, r, b, and g denote respectively $N$, $r$, $\beta$
% and $\gamma$, and the parameter F is a vector containing the values
% $f(\tilde s_{n+1/2}$
function I=singularintegral(N,r,b,g,F)
% nrows is the number of rows of tildeK1, tildeL1 tildeK2, tildeL2
nrows=max(2^ceil(log2(ceil(N/2)+N-1)),2);
% Compute $h_r=\pi/(2rN)$
hr=pi/(2*r*N);
% Define the functions $K_1(m,q)$, $K_2(m,q)$, $L_1(m,q)$ and $L_2(m,q)$ 
K1=@(m,q)(sin(hr*(q+1/2+2*r*m))./(hr*(q+1/2+2*r*m))).^b...
    .*((q+1+2*r*m).^(b+1)-(q+2*r*m).^(b+1)).*F(q+2*r*m+1);
L1=@(m,q)(sin(hr*(q+1/2-r-2*r*m))./(hr*(q+1/2-r-2*r*m))).^g...
    .*(sign(q+1-r-2*r*m).*abs(q+1-r-2*r*m).^(g+1)...
        -sign(q-r-2*r*m).*abs(q-r-2*r*m).^(g+1));
K2=@(m,q)(sin(hr*(r*N+q+1/2+2*r*m))./(hr*(r*N-q-1/2-2*r*m))).^b...
    .*((r*N-q-2*r*m).^(b+1)-(r*N-q-1-2*r*m).^(b+1)).*F(r*N+q+2*r*m+1);
L2=@(m,q)(sin(hr*(r*N+q+1/2-r-2*r*m))./(hr*(r*N+q+1/2-r-2*r*m))).^g...
    .*(sign(r*N+q+1-r-2*r*m).*abs(r*N+q+1-r-2*r*m).^(g+1)...
        -sign(r*N+q-r-2*r*m).*abs(r*N+q-r-2*r*m).^(g+1));
% Define the matrices $\tilde{\mathbf K}_1$, $\tilde{\mathbf K}_2$,
% $\tilde{\mathbf L}_1$, and $\tilde{\mathbf L}_2$
tildeK1=zeros(nrows,2*r);
tildeL1=zeros(nrows,2*r);
tildeK2=zeros(nrows,2*r);
tildeL2=zeros(nrows,2*r);
for q=0:r-1
    tildeK1(1:ceil(N/2),q+1)=K1((0:ceil(N/2)-1)',q);
    tildeK2(1:ceil(N/2),q+1)=K2((0:ceil(N/2)-1)',q);
    tildeL1([1:N end-ceil(N/2)+2:end],q+1)...
        =L1([0:N-1 -ceil(N/2)+1:-1]',q);
    tildeL2([1:N end-ceil(N/2)+2:end],q+1)...
        =L2([0:N-1 -ceil(N/2)+1:-1]',q);
end
for q=r:2*r-1
    tildeK1(1:floor(N/2),q+1)=K1((0:floor(N/2)-1)',q);
    tildeK2(1:floor(N/2),q+1)=K2((0:floor(N/2)-1)',q);
    tildeL1([1:N end-floor(N/2)+2:end],q+1)...
        =L1([0:N-1 -floor(N/2)+1:-1]',q);
    tildeL2([1:N end-floor(N/2)+2:end],q+1)...
        =L2([0:N-1 -floor(N/2)+1:-1]',q);
end
% Apply the fast convolution, to get the vector $\mathbf A$
A=ifft(sum(fft(tildeK1).*fft(tildeL1)+fft(tildeK2).*fft(tildeL2),2));
% Return the vector $\mathbf I$ containing the numerical approximation
% of the values $I(s_j)$
I=(hr^(b+g+1)/((b+1)*(g+1)))*A(1:N);