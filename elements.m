function v = elements(r,p,d)
Cd = r^2;
step = 0;
for k = 1:length(d)
    v(step+1:step+p(k)*Cd,1) = d(k)*ones(p(k)*Cd,1);
    step = step + p(k)*Cd;
end
end