function E = ImageEnergy(I,beta,nlp_in,nlp_out)
% Inputs:
% I = image being segmented
% beta = current contour
% nlp_in = negative log-likelihood of interior pixel values
% nlp_out = negative log-likelihood of exterior pixel values
% large = 1 if large image, 0 if small image (set to 0 unless image is very
% high-resolution)

% Output:
% E = image energy term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
N = length(nlp_in);

% If test image is large, do below in separate pieces to ensure no memory
% issues
if size(I,1)*size(I,2)>2e6
    large = 1;
else
    large = 0;
end

% Identify interior/exterior regions based on curve
mask = roipoly(I,beta(2,:),beta(1,:));
vmaski = mask(:);
vimg = I(:);
vmasko = logical(abs(vmaski-1));

% Add up negative log densities for interior/exterior
I_in = vimg(vmaski);
I_out = vimg(vmasko);

if large==1
    split = 4;         % to ensure enough memory to store below matrix
    li = length(I_in);
    lo = length(I_out);
    ssi = ceil(li/split);
    sso = ceil(lo/split);
    
    idx_i = []; idx_o = [];
    for i=1:(split-1)
        [~,tmp_i] = min(abs(repmat(I_in((1+(i-1)*ssi):i*ssi),1,N)-linspace(0,1,N)),[],2);
        [~,tmp_o] = min(abs(repmat(I_out((1+(i-1)*sso):i*sso),1,N)-linspace(0,1,N)),[],2);
        idx_i = [idx_i;tmp_i];
        idx_o = [idx_o;tmp_o];
    end
    [~,tmp_i] = min(abs(repmat(I_in(((split-1)*ssi+1):end),1,N)-linspace(0,1,N)),[],2);
    [~,tmp_o] = min(abs(repmat(I_out(((split-1)*sso+1):end),1,N)-linspace(0,1,N)),[],2);
    idx_i = [idx_i;tmp_i];
    idx_o = [idx_o;tmp_o];
else
    [~,idx_i] = min(abs(repmat(I_in,1,N)-linspace(0,1,N)),[],2);
    [~,idx_o] = min(abs(repmat(I_out,1,N)-linspace(0,1,N)),[],2);
end

E = sum(nlp_in(idx_i))+sum(nlp_out(idx_o));