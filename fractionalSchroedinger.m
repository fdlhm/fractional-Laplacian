% fractionalSchroedinger.m
% If you use it, please cite the corresponding paper:
% Jorge Cayama, Francisco de la Hoz, Carlota Mar\'{\i}a Cuesta, Carlos Javier Garc\'{\i}a-Cervera,
% A fast convolution method for the fractional Laplacian in $\mathbb{R}$, (2025).
%
% This program simulates the fractional cubic nonlinear Schr\"odinger equation in the focusing case
% by means of an explicit 4th order Runge-Kutta scheme.
clear
tic
L=200; % $L$
N=4096; % $N$
% sj is the vector containing the nodes $s_j$
sj=pi*((0:N-1)'+.5)/N;
% xj is the vector containing the nodes $x_j=\cot(s_j)$
xj=L*cot(sj);
% p0 is the vector containing the initial data $\psi_0(x)$
p0=exp(-xj.^2);
a=1.99; % $\alpha$
% rr is the vector containing all the values of $r$ to be considered
rr=2.^(0:6);
dt=1e-2; % $\Delta t$
mmax=1e4; % Number of iterations of the explicit RK4 scheme
mplot=1000; % Number of time instants to be stored
cplot=mmax/mplot; % cplot indicates how often a time instant will be stored
tall=dt*(0:mmax); % All the time instants
MM=zeros(mmax+1,length(rr)); % Energy of $\psi(x,t)$ for all $t$ and $r$
M=L*pi*sum(abs(p0).^2./sin(sj).^2)/N; % Energy of $\psi_0(x)$
MM(1,:)=M; % Store the energy of $\psi_0(x)$
pp=zeros(N,mplot+1); % Stored numerical approximations of $\psi(x,t)$
tt=zeros(mplot+1,1); % Stored time instants
for nr=1:length(rr)
    r=rr(nr); % Choose the corresponding $r$
    p=p0; % Initialize p with $\psi_0(x)$
    c=1; % Initialize the counter that indicates the array position of the results to be stored
    pp(:,c)=p0; % Store $\psi_0(x)$
    tt(c)=0; % Store the initial time
    c=c+1; % Update the counter
    for m=1:mmax
        k1u=-.5i*fractionalLaplacianfunction(p,r,a,L)+1i*abs(p).^2.*p;
        uaux=p+.5*dt*k1u;
        k2u=-.5i*fractionalLaplacianfunction(uaux,r,a,L)+1i*abs(uaux).^2.*uaux;
        uaux=p+.5*dt*k2u;
        k3u=-.5i*fractionalLaplacianfunction(uaux,r,a,L)+1i*abs(uaux).^2.*uaux;
        uaux=p+dt*k3u;
        k4u=-.5i*fractionalLaplacianfunction(uaux,r,a,L)+1i*abs(uaux).^2.*uaux;
        p=p+(dt/6)*(k1u+2*k2u+2*k3u+k4u); % Update $\psi$
        M=L*pi*sum(abs(p).^2 ./sin(sj).^2)/N; % Compute the energy
        MM(m+1,nr)=M; % Store the energy of $\psi$
        if M>2*MM(1,nr) % A simple stability criterion
            toc
            error('Stability error!')
        end
        if mod(m,cplot) == 0
            pp(:,c)=p; % Store $\psi$
            tt(c)=m*dt; % Store the time
            c=c+1; % Update the counter
        end
    end
    toc % Elapsed time
end
% Plot of the energy $M$
figure
h=plot(tall,log10(abs(MM-sqrt(pi/2))),'LineWidth',2);
set(gca,'YMinorGrid','off')
set(gca,'YMinorTick','off')
xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$\log_{10}(|M(t)-\sqrt{\pi/2}|)$','Interpreter','latex','FontSize',16)
title('$\alpha=1.99$, $L=200$, $N=4096$, $\Delta t=0.01$, $r\in\{2^0,2^1,\ldots,2^6\}$','Interpreter','latex','FontSize',16)
grid on
leg=cell(1,length(h));
for m=1:length(h)
    leg{m}=['$r=2^{',num2str(m-1),'}$'];
end
g=legend(leg,'Interpreter','latex','FontSize',16,'Location','southeast');
% Plot of the real part of $\psi(x,t)$
figure
mesh(tt,xj(1210:end-1209),real(pp(1210:end-1209,:)),sqrt(abs(pp(1210:end-1209,:))))
view([-45,45])
xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$x$','Interpreter','latex','FontSize',16)
zlabel('$\Re(\psi(x,t))$','Interpreter','latex','FontSize',16)
% Plot of the imaginary part of $\psi(x,t)$
figure
mesh(tt,xj(1210:end-1209),imag(pp(1210:end-1209,:)),sqrt(abs(pp(1210:end-1209,:))))
view([-45,45])
xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$x$','Interpreter','latex','FontSize',16)
zlabel('$\Im(\psi(x,t))$','Interpreter','latex','FontSize',16)
% Plot of the absolute value of $\psi(x,t)$
figure
mesh(tt,xj(1210:end-1209),abs(pp(1210:end-1209,:)),sqrt(abs(pp(1210:end-1209,:))))
view([-45,45])
xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$x$','Interpreter','latex','FontSize',16)
zlabel('$|\psi(x,t)|$','Interpreter','latex','FontSize',16)