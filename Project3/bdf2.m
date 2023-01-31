function [t,y,z] = bdf2(T,N,k,y0,z0,tol)
% T = terminal time
% N = number of time steps
% k = [R,C,w,A] parameter values
% y0 = initial conditions for [e2] as a column vector
% z0 = initial conditions for [e1;Iv] as a column vector
% tol = error tolerance for solver function

B=y0;
gam = 1;
y = y0;
z = z0;
dt = T/N;

tinput = dt;
[ytemp,ztemp] = daesolver(k,tinput,gam,dt,B,y(end),z(:,end),tol);
y = [y ytemp];
z = [z ztemp];
%B = [y(end);z(:,end)];
for n=2:N
    tinput = n*dt;
	% update y using BDF 2 method and your solver function to solve the implicit system.
    ytemp = (4/3)*y(:,end)-(1/3)*y(:,end-1);
    ztemp = (4/3)*z(:,end)-(1/3)*z(:,end-1);
    B = ytemp;
    gam = 2/3;
    [ytemp,ztemp]=daesolver(k,tinput,gam,dt,B,y(end),z(:,end),tol);
    y = [y ytemp];
    z = [z ztemp];
end
t = (0:N)*dt;
figure(1)
plot(t,y,'b-',t,z(1,:),'k-',t,z(2,:),'r-');
legend('e2 in Volts','e1 in Volts','Iv in Amperes');
xlabel('t (Seconds)')
ylabel('unit')
end