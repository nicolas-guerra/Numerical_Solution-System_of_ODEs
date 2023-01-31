%1b
x_0 = 1;
y_0 = 0;
gamma = 1;
T = 4*pi^2*(x_0^2+y_0^2)/gamma;
N = [50 100 200 400];

error = [];
figure(1)
for i = 1:length(N)
    dt = T/N(i);
    y = [y_0];
    x = [x_0];
    for j = 1:N(i)
        fy = x(end)/(2*pi*(x(end)^2+y(end)^2));
        fx = -y(end)/(2*pi*(x(end)^2+y(end)^2));        
        y = [y y(end)+dt*fy];
        x = [x x(end)+dt*fx];
    end
    error = [error sqrt((1-x(end))^2+(0-y(end))^2)];
    plot(x,y)
    hold on
end
%EXACT
th = 0:pi/50:2*pi;
xunit = 1 * cos(th);
yunit = 1 * sin(th);
plot(xunit, yunit);
title('Vortex Project 1(b)','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
legend('N = 50','N = 100','N = 200','N = 400','Exact','interpreter','latex')
hold off
disp(error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1c
error_rk = [];
figure(2)
for i = 1:length(N)
    [xrk,yrk,prk,qrk]=rk_nicolasguerra(N(i),T,[0],[0],[1],[1],[0]);
    error_rk = [error_rk sqrt((1-prk(end))^2+(0-qrk(end))^2)];
    plot(prk,qrk)
    hold on
end
%EXACT
th = 0:pi/50:2*pi;
xunit = 1 * cos(th);
yunit = 1 * sin(th);
plot(xunit, yunit);
title('Vortex Project 1(c)','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
legend('N = 50','N = 100','N = 200','N = 400','Exact','interpreter','latex')
hold off

disp(error_rk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2a
x = [0;.4;-0.2;-.2];
y = [0;0;.3464;-.3464];
g = [-1;-.667;-.667;-.667];
p = [1];
q = [0];
T = 20;
N = 400;

[x_out,y_out,p_out,q_out]=rk_nicolasguerra(N,T,x,y,g,p,q);
figure(3)
plot(p_out,q_out,'--')
hold on
for i = 1:size(x_out,1)
    plot(x_out(i,:),y_out(i,:))
end
title('Vortex Project 2(a)','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
legend('Particle 1','Vortex 1','Vortex 2','Vortex 3','Vortex 4','interpreter','latex')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2b
x = [0;1];
y = [0;0];
g = [1;-1];
p = [-.5;0;.5;1;1.5];
q = [.2;.2;.2;.2;.2];
T = 20;
N = 400;

[x_out,y_out,p_out,q_out]=rk_nicolasguerra(N,T,x,y,g,p,q);
figure(4)
for i = 1:size(p_out,1)
    plot(p_out(i,:),q_out(i,:),'--')
    hold on
end
for i = 1:size(x_out,1)
    plot(x_out(i,:),y_out(i,:))
end
title('Vortex Project 2(b)','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
legend('Particle 1','Particle 2','Particle 3','Particle 4','Particle 5','Vortex 1','Vortex 2','interpreter','latex')
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2c
x = [-.5;.5;0;-.5;.5;0];
y = [0;0;0;1;1;1];
g = [5;-5;.2;-5;5;-.2];
p = [0];
q = [.5];
T = 20;
N = 2000;

[x_out,y_out,p_out,q_out]=rk_nicolasguerra(N,T,x,y,g,p,q);
figure(5)
plot(p_out,q_out,'--')
hold on
for i = 1:size(x_out,1)
    plot(x_out(i,:),y_out(i,:))
end
title('Vortex Project 2(c)','interpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
legend('Particle 1','Vortex 1','Vortex 2','Vortex 3','Vortex 4','Vortex 5','Vortex 6','interpreter','latex')
hold off

