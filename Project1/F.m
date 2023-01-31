function [Fx,Fy,Fp,Fq] = F(x,y,g0,p,q)

Nv = max(size(x));
Np = max(size(p));
Fx = zeros(Nv,1);
Fy = zeros(Nv,1);
Fp = zeros(Np,1);
Fq = zeros(Np,1);

for j=1:Nv
	Fx(j) = 0;
	Fy(j) = 0;
	for k=1:Nv
		if j ~= k
			Fx(j) = Fx(j)-g0(k)/2/pi*(y(j)-y(k))/((x(j)-x(k))^2+(y(j)-y(k))^2);
			Fy(j) = Fy(j)+g0(k)/2/pi*(x(j)-x(k))/((x(j)-x(k))^2+(y(j)-y(k))^2);
		end
	end
end		
for j=1:Np
	Fp(j) = 0;
	Fq(j) = 0;
	for k=1:Nv
		Fp(j) = Fp(j)-g0(k)/2/pi*(q(j)-y(k))/((p(j)-x(k))^2+(q(j)-y(k))^2);
		Fq(j) = Fq(j)+g0(k)/2/pi*(p(j)-x(k))/((p(j)-x(k))^2+(q(j)-y(k))^2);
	end
end		

