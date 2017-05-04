function Wm = discrete_SVD_bb(y,B,resolution)
%   Input
%            y,             raw time domain samples per receiving antenna in one sub-frame, A-by-N,
%                               A=2,4,8,16,...  
%            B,             number of beams selected
%            resolution,    number of bits of digital phase shifters
%   Output
%            Wm,            beamforming weight, A-by-B

y=reshape(y,[size(y,1) length(y(:))./size(y,1)]);
[A,N]=size(y);
Rt=y*y'/N;
R = Rt;
dim = B;
M = size(Rt,1);
max_set = 1e6;
weight_d = exp(1j*[1:resolution]*2*pi/resolution);
Wm = zeros(M,dim);
Wm_rounded = zeros(M,dim);
for dim_index = 1:dim
    Q = eye(M) - Wm*Wm'/M;
    R = Q'*R*Q; 
    [w0,sigma0,~] = svd(R);
    w0 = w0(:,1);
    w0 = exp(1j*round(angle(w0)/(2*pi/resolution))*2*pi/resolution);
    w = w0;
    Wm_rounded(:,dim_index) = w;
    b = abs(w0'*R*w0/(w0'*w0));
    G = zeros(max_set,M);
    p = zeros(max_set,1);
    u = zeros(max_set,1);
    u(1) = sigma0(1,1);    
    G(1,:) = ones(1,M).*inf;
    G(1,1) = 1;
    p(1) = 1;
    set_index = 0;
    while(1)     
        [~,p_m] =  max(u);
%         if rand(1,1)>0.5;
%             [~,p_m] =  max(u);
% %             display('maxu')
%         else
%             [~,p_m] =  max(p);
% %             display('maxp')
%         end
        pt = p(p_m);
        gt = G(p_m,:);
        G(p_m,:) = zeros(1,M);
        p(p_m) = 0;
        u(p_m) = 0;
        gmax = find(abs(gt)>10);
        gmax = gmax(1);
        Gt = repmat(gt,[resolution,1]);
        Gt(:,gmax) = weight_d;
        pt = pt+1;    
        ut = zeros(resolution,1);
        eliminate_r = zeros(resolution,1);
        for i=1:resolution
            di = Gt(i,1:pt)';                
            [fval,ut(i),w_p] = partial_svd(R,di,B);
            if fval<b
                eliminate_r(i) = 1;
            else   
                if pt==M
                    eliminate_r(i) = 1;  
                    b = fval;
                    w = di;                      
                else
                    w_t = [di;w_p];
                    w_t = exp(1j*round(angle(w_t)/(2*pi/resolution))*2*pi/resolution);
                    b_t = abs(w_t'*R*w_t/(w_t'*w_t));
                    if (b_t>b)
                        b = b_t;
%                         display(b);
                        w = w_t;
                    end 
                end
            end
                 
        end
        Gt(logical(eliminate_r),:) = [];
        ut(logical(eliminate_r)) = [];
        rows = size(Gt,1);
        if rows~=0
            G((set_index+1):(set_index+rows),:) = Gt;
            u((set_index+1):(set_index+rows)) = ut;
            p((set_index+1):(set_index+rows)) = pt;
            set_index = set_index+rows;
        end
%         if(size(G,1)==0)
%             display(sprintf('optimization terminated for dim %d',dim_index));
%             break;
%         end
        umax = max(u);
        epi = 0.1*b;    
        if(umax-b<epi)
%             display(sprintf('optimization terminated (condition 1) for dim %d',dim_index));
            break;
        end 
    end
    Wm(:,dim_index) = w;
end
Wm = conj(Wm);

end

