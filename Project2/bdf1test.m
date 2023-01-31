function [t,y] = bdf1test(T,N,k,y0,tol)
% T = terminal time
% N = number of time steps
% k = [kd,ki,kt] parameters
% y0 = initial conditions for [I,M,R,Pdot,P] as a column vector

    gam = 1;
    y = y0;
    dt = T/N;
    C = y0;
    for n=1:N
    	% update y using BDF 1 method and your solver function to solve the implicit system.
        y = [y solver(k,gam,dt,C,y(:,end),tol)];
        C = y(:,end);
    end
    t=(0:N)*dt;
    


end
