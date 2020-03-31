function [z0,z,Relerrs] = solver(y,x, Params,A,At)
%% initialization
npower_iter = Params.npower_iter;

z0 = randn(Params.n1,Params.n2);
z0 = z0/norm(z0,'fro');

n  = Params.n1*Params.n2;
r  = Params.r;
m  = Params.m;
p  = Params.p;

normest = sqrt(sum(y(:).^2)/(n*r));
%     normest = sqrt(sum(y(:).^2)/numel(y(:)));

[~,B] = sort(y(:).^2,'descend');
ytr = zeros(size(y));
ytr(B(1:ceil((n*r)/p))) = 1;

for tt = 1:npower_iter                                                % Truncated power iterations
    z0 = At(ytr .* A(z0))./((n*r)*ceil((n*r)/p));
    z0 = z0/norm(z0,'fro');
end

z        = normest * z0;
Relerrs  = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro'); % Initial rel. error

u = Params.u0;

%% main loop
sgd    = Params.m;
batch  = Params.B;
gamma  = Params.y;
gamma1 = Params.y1;
mu = Params.mu;

for t = 1: Params.T
    for i=1:batch:sgd-batch+1
        
        ysub = zeros(size(y));
        ysub(i:i+batch-1) = 1;
        ayz  = A(z);
        yz   = sqrt(abs(ayz).^2+u^2);
        grad = At(ysub.*(ayz-y.*ayz./yz))/sgd;
        
        z = z - mu*grad;
        
        if norm(grad,'fro') < gamma*u
            u = gamma1*u;
        end
    end
    
    Relerrs(t + 1) = norm(x - exp(-1i*angle(trace(x'*z))) * z, 'fro')/norm(x,'fro');
    
end

end