function [t,sip] = rk_sip(sip0,param,T,N)
% sip0: [S0, I0, P0]
% param: [r, k, lambda, mu, m, a, theta, delta]
% T: terminal time
% N: number of time steps

dt = T/N;
sip = sip0;
r = param(1);
k = param(2);
lambda = param(3);
mu = param(4);
m = param(5);
a = param(6);
theta = param(7);
delta = param(8);

for step=1:(T/dt)
    
    [dS, dI, dP] = F_SIP(sip(1,step), sip(2,step), sip(3,step),param);
    K1 = dt*[dS; dI; dP];
    
    [dS,dI,dP] = F_SIP(sip(1,step) + K1(1)/2, sip(2,step)+K1(2)/2, sip(3,step)+ K1(3)/2,param);
    K2 = dt*[dS; dI; dP];
    
    [dS,dI,dP] = F_SIP(sip(1,step) + K2(1)/2, sip(2,step)+K2(2)/2, sip(3,step)+ K2(3)/2,param);
    K3 = dt*[dS; dI; dP];

    [dS,dI,dP] = F_SIP(sip(1,step) + K3(1), sip(2,step)+K3(2), sip(3,step)+ K3(3),param);
    K4 = dt*[dS; dI; dP];

    sip(:,step+1) = sip(:,step) + (K1+2*K2+2*K3+K4)/6.;
end
t = (0:N)*dt;
S = sip(1,:); I=sip(2,:); P=sip(3,:);
% include your solutions for the critical point here!
Sc = k-(delta*a*(k*lambda+r)/((theta-delta)*r));
Ic = delta*a/(theta-delta);
Pc = (a*theta*(r*(theta-delta)*(k*lambda-mu)-delta*a*lambda*(k*lambda+r)))/(r*m*(theta-delta)^2);


subplot(2,1,1);
plot(t,S,'r-',t,I,'g-',t,P,'b-',[t(1),t(end)],[Sc,Sc],'r:',[t(1),t(end)],[Ic,Ic],'g:',[t(1),t(end)],[Pc,Pc],'b:');
xlabel('Time')
ylabel('Population')
legend('S','I','P');
subplot(2,1,2)
plot3(S,I,P,'k-',Sc,Ic,Pc,'bo');
xlabel('S'); ylabel('I'); zlabel('P');
drawnow
end