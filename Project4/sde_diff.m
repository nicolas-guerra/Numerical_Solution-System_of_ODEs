function [t,S,I,P] = sde_diff(param,T,dt0,sip)
% Initial values sip: [S0, I0, P0]
% param: [r, k, lambda, mu, m, a, theta, delta]
% T: terminal time
% N: number of time steps

% intialize everything
sip0=sip;
S = sip0(1);
I = sip0(2);
P = sip0(3);
r = param(1);
k = param(2);
lambda = param(3);
mu = param(4);
m = param(5);
a = param(6);
theta = param(7);
delta = param(8);
t = 0;
n=1;

% take time-steps using Euler-Maruyama
while t (end) < T
    dt = min( [T-t(end), dt0 ]);
    dW = sqrt(dt)*randn(6, 1);
    % Compute S {n+1}, I {n+1}, P {n+1} using step size dt
    ac = [r*S(n)*(1-(S(n)+I(n))/k)-lambda*S(n)*I(n);...
      lambda*S(n)*I(n)-mu*I(n)-(m*I(n)*P(n))/(I(n)+a);...
      theta*I(n)*P(n)/(I(n)+a)-delta*P(n)];
  
    b = [sqrt(r*S(n))  -sqrt((r/k)*S(n)*(S(n)+I(n))) -sqrt(lambda*S(n)*I(n)) 0 0 0;...
         0 0 sqrt(lambda*S(n)*I(n)) -sqrt(m*P(n)*I(n)/(I(n)+a)+mu*I(n)) 0 0;...
         0 0 0 0 sqrt(theta*I(n)*P(n)/(I(n)+a))  -sqrt(delta*P(n))];
     
    S(n+1) = S(n)+ac(1)*dt+b(1,:)*dW;
    I(n+1) = I(n)+ac(2)*dt+b(2,:)*dW;
    P(n+1) = P(n)+ac(3)*dt+b(3,:)*dW;

        while min([S(n+1), I(n+1), P(n+1)]) < 0
        dt = dt/2 ;
        dW = dW/sqrt (2) ;

        % Recompute S {n+1}, I {n+1}, P {n+1} using new step size
        S(n+1) = S(n)+ac(1)*dt+b(1,:)*dW;
        I(n+1) = I(n)+ac(2)*dt+b(2,:)*dW;
        P(n+1) = P(n)+ac(3)*dt+b(3,:)*dW;
        end
    t = [t , t(end)+dt] ;
    n = n+1;

end

end