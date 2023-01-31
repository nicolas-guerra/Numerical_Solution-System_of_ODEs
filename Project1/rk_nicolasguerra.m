function [x_out,y_out,p_out,q_out] = rk_nicolasguerra(N, T, x, y, g, p, q)
    x_out = x;
    y_out = y;
    p_out = p;
    q_out = q;
    dt = T/N;
    for i = 1:N
        [Fx,Fy,Fp,Fq] = F(x_out(:,i),y_out(:,i),g,p_out(:,i),q_out(:,i));
        k1x = Fx;
        k1y = Fy;
        k1p = Fp;
        k1q = Fq;
        
        [Fx2,Fy2,Fp2,Fq2] = F(x_out(:,i)+(dt/2)*k1x,y_out(:,i)+(dt/2)*k1y,g,p_out(:,i)+(dt/2)*k1p,q_out(:,i)+(dt/2)*k1q);
        k2x = Fx2;
        k2y = Fy2;
        k2p = Fp2;
        k2q = Fq2;
        
        [Fx3,Fy3,Fp3,Fq3] = F(x_out(:,i)+(dt/2)*k2x,y_out(:,i)+(dt/2)*k2y,g,p_out(:,i)+(dt/2)*k2p,q_out(:,i)+(dt/2)*k2q);
        k3x = Fx3;
        k3y = Fy3;
        k3p = Fp3;
        k3q = Fq3;
        
        [Fx4,Fy4,Fp4,Fq4] = F(x_out(:,i)+dt*k3x,y_out(:,i)+dt*k3y,g,p_out(:,i)+dt*k3p,q_out(:,i)+dt*k3q);
        k4x = Fx4;
        k4y = Fy4;
        k4p = Fp4;
        k4q = Fq4;
        
        x_out(:,i+1) = x_out(:,i) + (dt/6)*(k1x+2*k2x+2*k3x+k4x);
        y_out(:,i+1) = y_out(:,i) + (dt/6)*(k1y+2*k2y+2*k3y+k4y);
        p_out(:,i+1) = p_out(:,i) + (dt/6)*(k1p+2*k2p+2*k3p+k4p);
        q_out(:,i+1) = q_out(:,i) + (dt/6)*(k1q+2*k2q+2*k3q+k4q);
    end
end