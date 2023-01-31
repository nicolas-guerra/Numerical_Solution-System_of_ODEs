t = 0:1/10:2*pi;
x = 1+2*cos(t);
y = 2*sin(t);
[theta, r] = cart2pol(x,y);

tol = 1e-14;
T = 10;

theta_bad = [];
r_bad = [];
for i=1:length(theta)
    x0 = [0,1,0,0];
    xT = [theta(i),r(i),0,0];
    
    [t,x,u]=spongebob(T,x0,xT);
    
    if norm(x(end,1:4)-xT,2)>tol
        theta_bad = [theta_bad theta(i)];
        r_bad = [r_bad r(i)];
    end
end

figure(1)
polarplot(theta,r, 'b')
hold on
polarplot(theta_bad,r_bad,'r*')
hold off