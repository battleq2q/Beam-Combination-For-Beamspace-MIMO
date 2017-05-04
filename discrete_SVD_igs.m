function Wm = discrete_SVD_igs(y,B,resolution)
% function [Wbf,eigenvalue]=svdBeamForming(y)
%            beamforming weights to combine uplink signals
%   Input
%            y,             raw time domain samples per receiving antenna in one sub-frame, A-by-N,
%                               A=2,4,8,16,...  
%            B,             number of beams selected
%            resolution,    number of bits of digital phase shifters
%            g,             number of weights serached each time
%   Output
%            Wbf,           beam forming weight, A-by-B
%            eigenvalue,    eigen values corresponding to each beam, B-by-1
%            power_ratio    power researved after BF
%

y=reshape(y,[size(y,1) length(y(:))./size(y,1)]);
[A,N]=size(y);
Rt=y*y'/N;
dim = B;
R = Rt;
M = size(Rt,1);
weight_d = exp(1j*[1:resolution]*2*pi/resolution);
Wm = zeros(M,dim);
for dim_index = 1:dim
    Q = eye(M) - Wm*Wm'/M;
    R = Q'*R*Q; 
    [w0,~,~] = svd(R);
    w0 = w0(:,1);
    w0 = exp(1j*round(angle(w0)/(2*pi/resolution))*2*pi/resolution);
    w = [];
    for i=1:M
        fval = 0;
        for j=1:resolution
            wt = [w;weight_d(j)];
            [fval(j),~,~] = partial_svd(R,wt,B);           
        end
        [~,iindex] = max(fval);
        w = [w;weight_d(iindex)];        
    end
    Wm(:,dim_index) = w;
end
Wm = conj(Wm);
end

