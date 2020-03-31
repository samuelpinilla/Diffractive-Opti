function D1 = spatial_distribution(n,p,d,r)

D1 = 2*ones(n); 
Cd = r^2;   
%% Guarantee the probability of the coding element 
v    = elements(r,p,d); 
ind  = randperm(Cd,Cd);
v1   = v(ind);
cel1 = reshape(v1,r,r);

%% Initial conditions
D1(1:r,1:r) = cel1;   
f0 = 1:r;  
c0 = 1:r;  
f1 = f0;
c1 = c0;
f2 = f1;
c2 = c1 + r;

%% spatial distribution 
for cel = 1:(n^2/Cd)-1   
    if mod(cel,n/r) == 0
        f1 = f0;
        c1 = c0;
        f2 = f1 + r;
        c2 = c1;
    elseif (cel~=1) && (mod(cel,n/r)==1)
        f0 = f2;
        c0 = c2;
        f1 = f0;
        f2 = f1;
        c1 = c0;
        c2 = c1 + r;
    elseif (cel~=1)
        f2 = f1;
        c1 = c1 + r;
        c2 = c1 + r;
    end
    
    aux = 0;
    for i1 = 1:length(f1)
        for j1 = 1:length(c1)
            
            dist = zeros(Cd,1);
            pos  = zeros(Cd,2);
            acum = 1;
            for i2 = 1:length(f2)
                for j2 = 1:length(c2)
                    
                    pos(acum,:) = [f2(i2),c2(j2)];
                    dist(acum)  = min(abs(f1(i1)-f2(i2)),abs(c1(j1)-c2(j2))); % set of distances 
                    acum = acum +1;
                end
            end
            
            if aux > 0
                for i = 1:size(pos,1)
                    if D1(pos(i,1),pos(i,2)) ~= 2
                        dist(i,:) = 0;
                    end
                end
            end
            
            [~,ind] = sort(dist,'descend');
            maximos = sum(dist == max(dist));
            ind = ind(1:maximos);  % positions where the maximum distances are found
            
            while 1
                temp = pos(ind(randperm(maximos,1)),:); % random position of the second cell
                if D1(temp(1),temp(2)) == 2
                    break
                end
            end
            D1(temp(1),temp(2)) = D1(f1(i1),c1(j1));
            aux = aux +1;
        end
    end
    aux = 0;
end
end


function v = elements(r,p,d)
    Cd = r^2;
    step = 0;
    for k = 1:length(d)
        v(step+1:step+p(k)*Cd,1) = d(k)*ones(p(k)*Cd,1);
        step = step + p(k)*Cd;
    end
end