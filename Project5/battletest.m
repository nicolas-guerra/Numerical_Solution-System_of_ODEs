function scr1 = battletest(a)

tf = 5+ceil(rand*5);
% theta0 = 2*pi*rand-pi;
% r0 = 1+5*rand;
% thetaf = theta0+2*pi*rand-pi;
% rf = 1+5*rand;
bc = 2*(2*rand(1,4)-1);
theta0 = atan2(bc(2),bc(1))
r0 = sqrt(bc(1)^2+bc(2)^2)
thetaf = atan2(bc(4),bc(3))
rf = sqrt(bc(3)^2+bc(4)^2)

fprintf(1,'\nChallenge:\nMove from (%f,%f) to (%f,%f) in time %d\n\n', ...
		r0*cos(theta0), r0*sin(theta0), rf*cos(thetaf), rf*sin(thetaf), tf);

fprintf('Running robot %s...\n',a);
tic;
[t1,x1,u1] = feval(a,tf,[theta0,r0,0,0],[thetaf,rf,0,0]);
time1 = toc;
% check energy compliance
nrgerror = max(max(abs([x1(:,7)./(1.5+5*x1(:,2).^2)-u1(:,1), x1(:,8)/5-u1(:,2)])));
if nrgerror > 1e-7
    [nrgerror,i] = max(max(abs([x1(:,7)./(1.5+5*x1(:,2).^2)-u1(:,1), x1(:,8)/5-u1(:,2)]')));
    [nrgerror,j] = max(max(abs([x1(:,7)./(1.5+5*x1(:,2).^2)-u1(:,1), x1(:,8)/5-u1(:,2)])));
    fprintf('There is an error in computing the energy vector (u1,u2)\n');
    if j==1
        fprintf('Compare: lambda3(%d)/(1.5+5*x2(%d)^2) = %e\n         u1(%d) = %e\n', ...
            i, i, x1(i,7)/(1.5+5*x1(i,2)^2), i, u1(i,1));
    else
        fprintf('Compare: lambda4(%d)/5 = %e\n         u2(%d) = %e\n', ...
            i, i, x1(i,8)/5, i, u1(i,2));
    end 
end
[ok,reason] = validate(t1,x1,u1,[theta0,r0,0,0],[thetaf,rf,0,0])
figure(1)
subplot(1,2,1);
nrg1 = plotrobot(t1,x1,u1);
title(a);
s = sprintf('Energy = %g',nrg1);
xlabel(s);
axis equal
drawnow
err1 = norm(x1(end,1:4)-[thetaf,rf,0,0]);
if isnan(err1)
    ax = axis;
    text((ax(1)+ax(2))/2,(ax(3)+ax(4))/2,'ERROR!');
end
subplot(1,2,2);
bar([1000*time1,nrg1,1000*err1]','stacked');
set(gca,'XTickLabel',{'time','energy','error'});
drawnow
scr1 = sum([1000*time1,nrg1,1000*err1]);

end

function energy = plotrobot(t,xx,u)

i = 1.5;
m = 5;
r1 = 1;
r2 = 1;

x = xx(:,2).*cos(xx(:,1));
y = xx(:,2).*sin(xx(:,1));
plot(x,y,'k-');
n = sqrt(x.^2+y.^2);

ax = x - u(:,1).*y./n/100;
ay = y + u(:,1).*x./n/100;

bx = x + u(:,2).*x./n/100;
by = y + u(:,2).*y./n/100;

hold on
for j=1:size(x,1)
    
    plot([x(j),ax(j)],[y(j),ay(j)],'b-');
    plot([x(j),bx(j)],[y(j),by(j)],'r-');
    
end
legend('extension force','torque force');
plot(0,0,'ko',0,0,'k+','MarkerSize',20);
hold off
axis equal

energy = trapz(t,u(:,1).^2+u(:,2).^2);
end