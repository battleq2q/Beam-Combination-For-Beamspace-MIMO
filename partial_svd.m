function [fval,fval_rd,w] = partial_svd( R,di,B )
d_dim = length(di);
M = size(R,1);
ml = M - d_dim;
if ml==0
    d = di'*di;
    r = di'*R*di;
    fval = r/d;
    fval_rd = fval;
    w = di(:);
    return;
end
di = di(:);
Ri = R(1:d_dim,1:d_dim);
Rj = R((d_dim+1):end,(d_dim+1):end);
Rji = R((d_dim+1):end,1:d_dim);
p = Rji * di;
[u,s,~] = svd(Rj);
pu = u'*p;
if abs(pu(1))<1e-8
    display('partial svd warning: corner case encountered');
end
d = di'*di;
r = di'*Ri*di;
lambda1 = s(1,1);
xl = lambda1;
xu = (d*lambda1+r+sqrt((d*lambda1-r)^2+4*d*p'*p))/2/d;
fval = abs(bi_sec(diag(s), pu, d, r, xl, xu));

w = (fval*eye(ml,ml)-Rj)\p;
x = [di;w];
sigma_e = 1/2^B*w'*w;
fval_rd = real((x'*R*x + sigma_e*trace(Rj)/ml)/(x'*x+sigma_e));  % rate-distortion theory correction
end

function x = bi_sec(lambda,pu,d,r,xl,xu)
x = (xu+xl)/2;
eps = 1e-2;
for i=1:10000
    if myf(x,lambda,pu,d,r)>0
        xl = x;
    else
        xu = x;
    end
    x = (xl+xu)/2;
    if(xu-xl < eps*lambda(1))
        break;
    end
end
if i==10000
    display('max iterations reached in bi-section');
end
end

function f = myf(x,lambda,pu,d,r)
s_vec = abs(pu).^2./(x-lambda);
f = sum(s_vec) - x*d+r;
end
