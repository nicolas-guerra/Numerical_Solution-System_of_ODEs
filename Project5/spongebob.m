function [t,x,u]=spongebob(T,x0,xT)
%Written entirely and originally by Nicolas Guerra

%guess s
s = [0,0,0,0];
y0 = [x0,s];
N=2000;
tol = 1e-14;
m=5;
I = 1.5;

[t_temp,y_temp]=rk(y0,T,N);
%calculate residuals
h = [y_temp(end,1)-xT(1);y_temp(end,2)-xT(2);...
    y_temp(end,3)-xT(3);y_temp(end,4)-xT(4)];

%check if residual is good enough
while norm(h,2)>tol
    %dh is the final values of the numerically solved z1...4
    dh = [y_temp(end,9)  y_temp(end,17) y_temp(end,25) y_temp(end,33);
          y_temp(end,10) y_temp(end,18) y_temp(end,26) y_temp(end,34);
          y_temp(end,11) y_temp(end,19) y_temp(end,27) y_temp(end,35);
          y_temp(end,12) y_temp(end,20) y_temp(end,28) y_temp(end,36)];
    if rcond(dh)<1e-15
        break
    end
    %update s
    s = s - (dh\h)';
    y0=[x0,s];
    [t_temp,y_temp]=rk(y0,T,N);
    h = [y_temp(end,1)-xT(1);y_temp(end,2)-xT(2);...
    y_temp(end,3)-xT(3);y_temp(end,4)-xT(4)];
end
t = t_temp;
x = y_temp(:,1:8);
u1 = x(:,7)./(I+m*x(:,2).^2);
u2 = x(:,8)./m;
u = [u1 u2];
    
end

function [t, y]= rk(y0,T,N)
%y is the current values of x1...4 and lambda1...4
%y = [x1,x2,x3,x4,lam1,am2,lam3,lam4]
dt = T/N;

x1=[y0(1)];
x2=[y0(2)];
x3=[y0(3)];
x4=[y0(4)];
lam1=[y0(5)];
lam2=[y0(6)];
lam3=[y0(7)];
lam4=[y0(8)];
z1 = [0,0,0,0,1,0,0,0];
z2 = [0,0,0,0,0,1,0,0];
z3 = [0,0,0,0,0,0,1,0];
z4 = [0,0,0,0,0,0,0,1];

y = [y0 z1 z2 z3 z4];
%y = [y0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for step=1:N
    k1 = F(x1(end),x2(end),x3(end),x4(end),...
        lam1(end),lam2(end),lam3(end),lam4(end),...
        z1',z2',z3',z4');
    k1 = dt*k1;
    k2 = F(x1(end)+k1(1)/2,x2(end)+k1(2)/2,x3(end)+k1(3)/2,...
        x4(end)+k1(4)/2,lam1(end)+k1(5)/2,lam2(end)+k1(6)/2,...
        lam3(end)+k1(7)/2,lam4(end)+k1(8)/2,...
        z1'+k1(9)/2,z2'+k1(10)/2,z3'+k1(11)/2,z4'+k1(12)/2);
    k2 = dt*k2;
    k3 = F(x1(end)+k2(1)/2,x2(end)+k2(2)/2,x3(end)+k2(3)/2,...
        x4(end)+k2(4)/2,lam1(end)+k2(5)/2,lam2(end)+k2(6)/2,...
        lam3(end)+k2(7)/2,lam4(end)+k2(8)/2,...
        z1'+k2(9)/2,z2'+k2(10)/2,z3'+k2(11)/2,z4'+k2(12)/2);
    k3 = dt*k3;
    k4 = F(x1(end)+k3(1),x2(end)+k3(2),x3(end)+k3(3),...
        x4(end)+k3(4),lam1(end)+k3(5),lam2(end)+k3(6),...
        lam3(end)+k3(7),lam4(end)+k3(8),...
        z1'+k3(9),z2'+k3(10),z3'+k3(11),z4'+k3(12));
    k4 = dt*k4;
    
    y = [y; y(end,:)+((k1+2*k2+2*k3+k4)/6)'];  
    x1=[x1 y(end,1)];
    x2=[x2 y(end,2)];
    x3=[x3 y(end,3)];
    x4=[x4 y(end,4)];
    lam1=[lam1 y(end,5)];
    lam2=[lam2 y(end,6)];
    lam3=[lam3 y(end,7)];
    lam4=[lam4 y(end,8)];
    z1 = y(end,9:16);
    z2 = y(end,17:24);
    z3 = y(end,25:32);
    z4 = y(end,33:40);
end
    t = (0:N)*dt;

end

function f=F(x1,x2,x3,x4,lam1,lam2,lam3,lam4,z1,z2,z3,z4)
m=5;
I = 1.5;

x1dot = x3;
x2dot = x4;
x3dot = (1/(I+m*x2^2))*(lam3/(I+m*x2^2)-2*m*x2*x3*x4);
x4dot = lam4/m^2+x2*x3^2;
lam1dot = 0;
lam2dot = (2*m/(I+m*x2^2))*(x3*x4+(x2/(I+m*x2^2))*(lam3/(I+m*x2^2)-2*m*x2*x3*x4))*lam3-x3^2*lam4;
lam3dot = -lam1+2*m*x2*x4*lam3/(I+m*x2^2)-2*x2*x3*lam4;
lam4dot = -lam2+2*m*x2*x3*lam3/(I+m*x2^2);

%jacobian of the original 8 ODEs. These will make the odes for z.
j = [0,0,1,0,0,0,0,0;
     0,0,0,1,0,0,0,0;
     0,(((-2*lam3*m*x2/(I+m*x2^2)^2)-2*m*x3*x4)/(I+m*x2^2))-(2*m*x2*(lam3/(I+m*x2^2)-2*m*x2*x3*x4)/(I+m*x2^2)^2), -2*m*x2*x4/(I+m*x2^2),-2*m*x2*x3/(I+m*x2^2),0,0,1/(I+m*x2^2)^2,0;
     0,x3^2,2*x2*x3,0,0,0,0,1/m^2;
     0,0,0,0,0,0,0,0;
     0,(2*lam3*m*(lam3*(I-5*m*x2^2)+2*m*x2*(3+m*x2^2*(-2*I+m*x2^2))*x3*x4))/(I+m*x2^2)^4,-2*lam4*x3-(2*lam3*m*(-I+m*x2^2)*x4)/(I+m*x2^2)^2,-2*lam3*m*(-I+m*x2^2)*x3/(I+m*x2^2)^2,0,0,-2*m*(-2*lam3*x2+(I+m^2*x2^4)*x3*x4)/(I+m*x2^2)^3,-x3^2;
     0,-2*lam4*x3-2*lam3*m*(-I+m*x2^2)*x4/(I+m*x2^2)^2,-2*lam4*x2,2*lam3*m*x2/(I+m*x2^2),-1,0,2*m*x2*x4/(I+m*x2^2),-2*x2*x3;
     0,-2*lam3*m*(-I+m*x2^2)*x3/(I+m*x2^2)^2,2*lam3*m*x2/(I+m*x2^2),0,0,-1,2*m*x2*x3/(I+m*x2^2),0];

dz1 = j*z1;
dz2 = j*z2;
dz3 = j*z3;
dz4 = j*z4;


f = [x1dot;x2dot;x3dot;x4dot;...
    lam1dot;lam2dot;lam3dot;lam4dot;...
    dz1;dz2;dz3;dz4];

end