function [t,y,z] = dae(T,idt,k,y0,z0,tol,itol)
% T = final time
% idt = initial time step
% k = [R,C,w,A] parameter values
% y0=initial value for y, and should be the previous time step
% z0=initial values for z, and should be the previous time step
% tol=global tolerance
% itol = tolerance for the Newton solver

y = y0;
z = z0;

dt = idt;
gam = 1;
tinput = dt;
B = y0;
[ytemp,ztemp] = daesolver(k,tinput,gam,dt,B,y(end),z(:,end),itol);
y = [y ytemp];
z = [z ztemp];

t = [0, dt];
while t(end)<T
    dtn = dt;
    dtnm1 = t(end)-t(end-1);
    
    tinput = t(end)+dtn;
    gam = 1;
    B = y(end);
    [y1,z1] = daesolver(k,tinput,gam,dtn,B,y(end),z(:,end),itol);
    gam = (dtn+dtnm1)/(2*dtn+dtnm1);
    B = (dtnm1+dtn)^2/dtnm1/(2*dtn+dtnm1)*y(end)-dtn^2/dtnm1/(2*dtn+dtnm1)*y(end-1);
    [y2,z2] = daesolver(k,tinput,gam,dtn,B,y(end),z(:,end),itol);
    
    err = norm([y1;z1]-[y2;z2],2);
    dtstar = 0.9*tol/T*dtn^2/err;
    dtn = min([dtstar, 2*dt, T-t(end)]);
    
    tinput = t(end)+dtn;
    gam = 1;
    B = y(end);
    [y1,z1] = daesolver(k,tinput,gam,dtn,B,y(end),z(:,end),itol);
    gam = (dtn+dtnm1)/(2*dtn+dtnm1);
    B = (dtnm1+dtn)^2/dtnm1/(2*dtn+dtnm1)*y(end)-dtn^2/dtnm1/(2*dtn+dtnm1)*y(end-1);
    [y2,z2] = daesolver(k,tinput,gam,dtn,B,y(end),z(:,end),itol);
    
    while norm([y1;z1]-[y2;z2],2)>dtn/T*tol
        err = norm([y1;z1]-[y2;z2],2);
        dtstar = 0.9*tol/T*dtn^2/err;
        dtn = min([dtstar, 2*dt, T-t(end)]);
        
        tinput = t(end)+dtn;
        gam = 1;
        B = y(end);
        [y1,z1] = daesolver(k,tinput,gam,dtn,B,y(end),z(:,end),itol);
        gam = (dtn+dtnm1)/(2*dtn+dtnm1);
        B = (dtnm1+dtn)^2/dtnm1/(2*dtn+dtnm1)*y(end)-dtn^2/dtnm1/(2*dtn+dtnm1)*y(end-1);
        [y2,z2] = daesolver(k,tinput,gam,dtn,B,y(end),z(:,end),itol);
    end
    dt = dtn;
    t = [t,t(end)+dtn];
    y = [y y1];
    z = [z z1];
end
figure(1)
subplot(2,1,1)
plot(t,y,'b-',t,z(1,:),'k-',t,z(2,:),'r-');
legend('e2 in Volts','e1 in Volts','Iv in Amperes');
xlabel('t (Seconds)')
ylabel('unit')
subplot(2,1,2)
plot(t(1:end-1),t(2:end)-t(1:end-1),'k-')
xlabel('t (Seconds)')
ylabel('dt')
end