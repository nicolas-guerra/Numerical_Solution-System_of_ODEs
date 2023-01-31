function [t,y] = bdf2test(T,N,k,y0,tol)
% y = bdf1test(T,dt,k,y0,tol)
% T = terminal time
% N = number of time steps
% k = [kd,ki,kt] parameters
% y0 = initial conditions for [I,M,R,Pdot,P] as a column vector
% tol = error tolerance for solver function
C=y0;
gam = 1;
y = y0;
dt = T/N;

y = [y solver(k,gam,dt,C,y(:,end),tol)];
C = y(:,end);

for n=2:N
	% update y using BDF 2 method and your solver function to solve the implicit system.
    y_temp = (4/3)*y(:,end)-(1/3)*y(:,end-1);
    C = y_temp;
    y = [y solver(k,2/3,dt,C,y(:,end),tol)];
end
t=(0:N)*dt;

end
