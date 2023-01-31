function J = JacF(k,dt)
% k = [R,C,w,A] parameter values
kr = k(1);
kc = k(2);

% dt is the delta time
J = [1+(dt/(kr*kc)) -dt/(kr*kc) 0;
     -1/kr         1/kr      -1;
      0            1          0];
end