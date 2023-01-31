function [dS, dI, dP] = F_SIP(S,I,P,param)
r = param(1);
k = param(2);
lambda = param(3);
mu = param(4);
m = param(5);
a = param(6);
theta = param(7);
delta = param(8);


dS = r*S*(1-(S+I)/k)-lambda*S*I;
dI = lambda*S*I-mu*I-m*I*P/(I+a);
dP = theta*I*P/(I+a) - delta*P;

end
