function [y,z] = daesolver(k,t,gam,dt,B,yn,zn,tol)

% k = [R,C,w,A] parameter values
% t is t_(n+1) when solving for y_(n+1)
% gam=gamma
% dt=time step size
% B=constant vector
% yn=initial value for y, and should be the previous time step
% zn=initial values for z, and should be the previous time step
% tol=tolerance for the error in the residual
kr = k(1);
kc = k(2);
kw = k(3);
ka = k(4);

J = JacF(k,dt);
%initialize z and compute initial residual r
p = yn; % yn = [e2]
q = zn; % zn = [e1; Iv]
r = [p-gam*dt*(zn(1)-yn)/(kr*kc)-B;
    -zn(2)+(zn(1)-yn)/kr;
    ka*sin(kw*t)+zn(1)];

%Build the matrix A = I-dt*gamma*J
pq = [p;q];
%A = eye(3)-dt*gam*J;
while norm(r,2)>tol
    %update the value of z. you may use the MATLAB backslash
    pq = [p;q] - J\r;
    p = pq(1);
    q = pq(2:3);
    %update the value of the residual

r = [p-gam*dt*(q(1)-p)/(kr*kc)-B;
    -q(2)+(q(1)-p)/kr;
    ka*sin(kw*t)+q(1)];
end
y = pq(1);
z = pq(2:3);
 
end
