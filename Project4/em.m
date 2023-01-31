function [t,S,I,P] = em(param,T,N,sip0)
% Initial values sip: [S0, I0, P0]
% param: [r, k, lambda, mu, m, a, theta, delta]
% T: terminal time
% N: number of time steps

% intialize everything
dt = T/N;
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

% take time-steps using Euler-Maruyama
%(do not correct for negative values)
for n = 1:N
    W = [randn;randn;randn;randn;randn;randn];
    
    ac = [r*S(n)*(1-(S(n)+I(n))/k)-lambda*S(n)*I(n);...
      lambda*S(n)*I(n)-mu*I(n)-(m*I(n)*P(n))/(I(n)+a);...
      theta*I(n)*P(n)/(I(n)+a)-delta*P(n)];
  
    b = [sqrt(r*S(n))  -sqrt((r/k)*S(n)*(S(n)+I(n))) -sqrt(lambda*S(n)*I(n)) 0 0 0;...
         0 0 sqrt(lambda*S(n)*I(n)) -sqrt(m*P(n)*I(n)/(I(n)+a)+mu*I(n)) 0 0;...
         0 0 0 0 sqrt(theta*I(n)*P(n)/(I(n)+a))  -sqrt(delta*P(n))];
    
    S(n+1) = S(n)+ac(1)*dt+(b(1,:)*W)*sqrt(dt);
    I(n+1) = I(n)+ac(2)*dt+(b(2,:)*W)*sqrt(dt);
    P(n+1) = P(n)+ac(3)*dt+(b(3,:)*W)*sqrt(dt);
end
t = (0:N)*dt;
end