function [ok,reason,err_score] = validate(t,x,u,x0,xT)

if nargin == 1
    x0 = [0,1,0,0];
    xT = [pi/2,1,0,0];
    [t,x,u] = feval(t,3,x0,xT);
end

I = 1.5;
m = 5;

% validate x has the correct dimensions
dimok = size(x,2) == 8 && size(x,1) > 2 && size(u,1) == size(x,1) ...
    && length(t) == size(x,1);

% validate x values against differential equation
if dimok
    odeerr = max(max(abs(G(t,x,x0,xT))));
else
    odeerr = 0;
end

% validate u vectors
if dimok
    nrgerr = max(max(abs([x(:,7)./(1.5+5*x(:,2).^2)-u(:,1), x(:,8)/5-u(:,2)])));
    if nrgerr > 1e-7
        [nrgerr,i] = max(max(abs([x(:,7)./(1.5+5*x(:,2).^2)-u(:,1), x(:,8)/5-u(:,2)]')));
        [nrgerr,j] = max(max(abs([x(:,7)./(1.5+5*x(:,2).^2)-u(:,1), x(:,8)/5-u(:,2)])));
        fprintf('There is an error in computing the energy vector (u1,u2)\n');
        if j==1
            fprintf('Compare: lambda3(%d)/(1.5+5*x2(%d)^2) = %e\n         u1(%d) = %e\n', ...
                i, i, x(i,7)/(1.5+5*x(i,2)^2), i, u(i,1));
        else
            fprintf('Compare: lambda4(%d)/5 = %e\n         u2(%d) = %e\n', ...
                i, i, x(i,8)/5, i, u(i,2));
        end 
    end
else
    nrgerr = 0;
end

err_score = odeerr + nrgerr;
ok = (odeerr < 1e-1) && (nrgerr < 1e-5) && dimok;
if ok == 0
    if dimok == 0
        if size(x,2) ~= 8 || size(x,1) <=2
            reason = sprintf('Dimensions of x are %dx%d\n',size(x,1),size(x,2));
        elseif size(u,1) ~= size(x,1)
            reason = sprintf('Dimensions of x are %dx%d, but the dimensions of u are %dx%d\n', ...
                size(x,1),size(x,2),size(u,1),size(u,2));
        else
            reason = sprintf('Dimensions of x are %dx%d, but t has length %d\n', ...
                size(x,1),size(x,2),length(t,1));
        end
    elseif odeerr >= 1e-1
        reason = sprintf('ODE relative error too large: %e\nNeeds to be < 1.0e-1, see Figure(1)\n',odeerr);
    else
        reason = sprintf('U vector error too large: %e\nCheck computation of u\n',nrgerr);
    end
    fprintf('Not OK, %s',reason);
else
    fprintf('Code checks OK, Error score is %e\n', err_score);
    reason = '';
end
end

function g = G(t,y,x0,xT)
[tt,yy] = ode45(@F,[t(1),t(end)],y(1,:));
g = [];
figure(1)
ls=['r-','g-','b-','y-','r--','g--','b--','y--'];
sym=['r:o','g:+','b:*','y:x','r:s','g:d','b:^','y:v'];
vrbl={'x_1','x_2','x_3','x_4','lambda_1','labmda_2','lambda_3','lambda_4'};
for i=1:8
    ycheck = spline(tt,yy(:,i),t);
    h = (ycheck'-y(:,i))/max(abs(ycheck));
    g = [g, h];
    figure(1)
    subplot(4,2,i)
    plot(t,y(:,i),'b--',tt,yy(:,i),'r-');
    legend(sprintf('your %s',vrbl{i}),sprintf('exact %s',vrbl{i}));
end
%g = [yy(1,1:4)-x0, yy(end,1:4)-xT];
end

function f = F(t,y)
m=5;
I=1.5;
d=I+m*y(2)^2;
f = [y(3); ...
    y(4); ...
    1/d*(y(7)/d-2*m*y(2)*y(3)*y(4)); ...
    y(8)/m^2+y(2)*y(3)^2; ...
    0; ...
    2*m/d*(y(3)*y(4)+y(2)/d*(y(7)/d-2*m*y(2)*y(3)*y(4)))*y(7)-y(3)^2*y(8); ...
    -y(5)+2*m/d*y(2)*y(4)*y(7)-2*y(2)*y(3)*y(8); ...
    -y(6)+2*m/d*y(2)*y(3)*y(7)];
end

