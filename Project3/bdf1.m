function [t,y,z] = bdf1(T,N,k,y0,z0,tol)
% T = terminal time
% N = number of time steps
% k = [R,C,w,A] parameter values
% y0 = initial conditions for [e2] as a column vector
% z0 = initial conditions for [e1;Iv] as a column vector
% tol = error tolerance for solver function

gam = 1;

y = y0;
z = z0;
B = y0;
dt = T/N;
for n=1:N
    t = n*dt;
    [ynew,znew] = daesolver(k,t,gam,dt,B,y(:,end),z(:,end),tol);
    y = [y ynew];
    z = [z znew];
    B = y(:,end);
end
t = (0:N)*dt;
figure(1)
plot(t,y,'b-',t,z(1,:),'k-',t,z(2,:),'r-');
legend('e2 in Volts','e1 in Volts','Iv in Amperes');
xlabel('t (Seconds)')
ylabel('unit')
end