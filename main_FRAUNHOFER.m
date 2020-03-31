clc;
clear;
close all;

%% parameters 
addpath('Designed_Apertures')
% designed
cod = {};
cod{1} = [1,0];                                % cod1 
cod{2} = [1,0,1i];                             % cod2
cod{3} = [1,-1,1i,-1i];                        % cod3

% probability of designed
proba = {};
proba{1} = [1/2,1/2];            
proba{2} = [1/3,1/3,1/3];      
proba{3} = [1/4,1/4,1/4,1/4];     

% random
cod{4} = [1,0];                                 % cod1 
cod{5} = [1,0,1i];                              % cod2
cod{6} = [1,-1,1i,-1i];                         % cod3

%% variables
l           = 12; % shots
n_test      = 50;
results_mat = zeros(size(cod,2),l);

for ii = 1:size(cod,2)
    if ii == 2 || ii == 5
        N       = 63;
        ss      = 3;
    else
        N       = 64;
        ss      = 4;
    end
    x       = randn(N) + 1i*randn(N);
    [n1,n2] = size(x);
    Params.n1  = n1; 
    Params.n2  = n2; 
    
    for L = 1:l
        if exist('Params')                == 0,  Params.n1          = n1;               end
        if isfield(Params, 'n2')          == 0,  Params.n2          = n2;               end
        if isfield(Params, 'T')           == 0,  Params.T           = 250;              end
        if isfield(Params, 'y1')          == 0,  Params.y1          = 0.6;              end
        if isfield(Params, 'p')           == 0,  Params.p           = 2;                end
        if isfield(Params, 'u0')          == 0,  Params.u0          = 30;               end
        if isfield(Params, 'y')           == 0,  Params.y           = 0.1;              end
        if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 800;              end
        if isfield(Params, 'mu')          == 0,  Params.mu          = 4.2;              end
        if isfield(Params, 'alpha')       == 0,  Params.alpha       = 0.5;              end
        if isfield(Params, 'B')           == 0,  Params.B           = (n1*n2)/2;        end
        
        Params.L    = L;
        m           = floor(n1*n2*L);
        Params.m    = m;
        
        prob = 0;
        for hh = 1:n_test
            if ii<=3
                Masks = codes(N,L,proba{ii},cod{ii},ss);
            else
                Masks   = zeros(n1,n2,L);
                for ll = 1:L, Masks(:,:,ll) = randsrc(n1,n2,cod{ii}); end
            end
            
            r            = mean(mean(sum(abs(Masks).^2,3)));
            Params.r = r;
            Params.Masks = Masks;
            Params.ii    = ii;
            
            display(Params)
                     
            % Make linear operators;
            A  = @(I) fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));
            At = @(I) sum(Masks .* ifft2(I), 3) * size(I,1) * size(I,2);

            %% Make signal and data (noiseless)
            y = abs(A(x));
            
            [z0,z,Relerrs] = solver(y,x,Params, A,At);
            
            %% results
            prob = prob + (Relerrs(end) < 1e-5);
            
        end
        prob = prob /n_test;
        results_mat(ii,L) = prob;
    end
end
save('result_fraunhofer_success','results_mat');
