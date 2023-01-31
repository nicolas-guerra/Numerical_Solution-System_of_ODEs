function z = solver(k,gam,dt,C,yn,tol)
% z = solver(k,gam,dt,C,yn,tol)
% solves z = gamma*dt*F(z) + C
% where
% k = [kd, ki, kt] parameter values
% gam = gamma 
% dt = time step size
% C = constant vector
% yn = initial value for z, and should be the previous time step
% tol = tolerance for the error in the residual.

J = JacF(k,yn);

% initialize z and compute initial residual r
z = yn;
r = z-gam*dt*F(k,z)-C;
% Build the matrix A = I - dt*gamma*J
A = eye(5)-dt*gam*J;
while norm(r,2) > tol
    % update the value of z. You may use the Matlab backslash operator to compute A^(-1)r
    z = z - A\r;
    % update the value of the residual
    r = z-gam*dt*F(k,z)-C;
end

end % end of your solver function


% subfunction for F (allows functions in the same file to have access to F)
function yp = F(k,y)
% k = [kd, ki, kt] parameter values
kd = k(1);
ki = k(2);
kt = k(3);
I = y(1);
R = y(2);
M = y(3);
Pdot = y(4);
P = y(5);
% return the value of F as a column vector of length 5
yp = [-kd*I;
    2*kd*I-ki*M*R-2*kt*R^2-kt*R*Pdot;
    -ki*M*R-ki*M*Pdot;
    ki*M*R-kt*R*Pdot-2*kt*Pdot^2;
    kt*Pdot*R+kt*Pdot^2];
end

% subfunction for JacF (allows functions in the same file to have access to JacF)
function J = JacF(k,y)
% k = [kd, ki, kt] parameter values
kd = k(1);
ki = k(2);
kt = k(3);
I = y(1);
R = y(2);
M = y(3);
Pdot = y(4);
P = y(5);
% return the Jacobian dF/dy as a 5x5 matrix
% evaluated at y
J = [-kd 0 0 0 0;
    2*kd -ki*M-4*kt*R-kt*Pdot -ki*R -kt*R 0;
    0 -ki*M -ki*R-ki*Pdot -ki*M 0;
    0 ki*M-kt*Pdot ki*R -kt*R-4*kt*Pdot 0;
    0 kt*Pdot 0 kt*R+2*kt*Pdot 0];
end
