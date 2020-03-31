%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design of coded apertures 

%%% Input: 
    % n : coded aperture size (square) 
    % L : # shots  
    % r : cell size (square), Cd = r^2 cardinality (# pixels in a cell)
    % d : random variable 
    % p : probability 
    
%%% Output: 
    % D : Coded aperture designed 
    
%%% Auxiliar functions: 
    % elements.m 
    % spatial_distribution.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = codes(n,L,p,d,r)

% Generates the first shot
D1 = spatial_distribution(n,p,d,r);
D  = cat(3,D1,2*ones(n,n,L-1));

% from second shot 
for l = 2:L
    daux = elements(r,p,d);
    i  = 1;
    j  = 1;
    i0 = 1;
    j0 = 1;
    
    while sum(sum(D(:,:,l)==2)) > 0
        dif  = daux(D1(i,j) ~= daux); 
        if isempty(dif)
            D(i,j,l) = daux(1);
            daux(1) = [];
        else
            temp = randperm(length(dif),1);
            %                 [~,temp] = max(dif(:,2));
            D(i,j,l) = dif(temp);
            elim     = find(daux == dif(temp));
            daux(elim(1)) = [];
        end
        
        % shot ready! 
        if i == n && j == n
            break
        end
        
        % cell iterations
        if mod(i,r) == 0 && mod(j,r) == 0 && j ==n
            daux = elements(r,p,d);
            i0 = i0 + r;
            i  = i0;
            j0 = 1;
            j  = j0;
        elseif mod(i,r) == 0 && mod(j,r) == 0
            daux = elements(r,p,d);
            j0 = j0 + r;
            j  = j0;
            i  = i0;
        elseif mod(j,r) == 0
            i = i + 1;
            j = j0;
        else
            j = j + 1;
        end
    end
    
    % If d is binary 
    if sum(d==0)== 1 && size(d,2) == 2 
        % If the shot-th is even 
        if (-1)^l == 1 
            D1 = spatial_distribution(n,p,d,r);
        else
            D1 = D(:,:,l);
        end
    else
        D1 = D(:,:,l);
    end
end
end
    

    