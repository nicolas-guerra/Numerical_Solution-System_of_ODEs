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

% Build the matrix A = I - dt*gamma*J

while norm(r,2) > tol
    % update the value of z. You may use the Matlab backslash operator to compute A^(-1)r

    % update the value of the residual

end

end % end of your solver function


% subfunction for F (allows functions in the same file to have access to F)
function yp = F(k,y)
% k = [kd, ki, kt] parameter values

% return the value of F as a column vector of length 5

end

% subfunction for JacF (allows functions in the same file to have access to JacF)
function J = JacF(k,y)
% k = [kd, ki, kt] parameter values

% return the Jacobian dF/dy as a 5x5 matrix
% evaluated at y

end
