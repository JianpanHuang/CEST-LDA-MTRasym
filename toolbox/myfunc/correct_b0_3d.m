function img_corr = correct_b0_3d(img, db0, offs, roi, method)
% FUNCTION:
%   To correct B0 for images (CEST)
% INPUT:
%   img: images used for deltaB0 map generation
%   db0: deltaB0 map
%   offs: frequency offsets
%   roi: ROI for deltaB0 map correction
%   method: spline, makima or pchip
% OUTPUT:
%   db0: images with B0 correction
% REFERENCE
%   Kim M, et al. Magnetic Resonance in Medicine 2009, 61(6): 1441-1450.
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
if nargin < 4
    roi = ones(size(img));
end 
if nargin < 5
    method = 'spline';
end

% correct B0
sz = size(img);
if length(sz)<4
    sz(4) = sz(3);
    sz(3) = 1;
    img = reshape(img,sz);
    roi = reshape(roi,sz);
end
img_corr = zeros(sz);
for s = 1:sz(3)
    for m = 1:sz(1)
        for n = 1:sz(2)
            if roi(m,n,s) == 1
                db0_pix = db0(m,n,s);                 
                z_spec = squeeze(img(m,n,s,:));
                z_spec_max = max(z_spec(:));
                [~, off_max_ind] = max(offs);
                z_spec_corr = interp1(offs-db0_pix, z_spec, offs, method, z_spec(off_max_ind));
                z_spec_corr = z_spec_corr/((max(z_spec_corr(:))/(z_spec_max+eps))+eps);
                % z_spec_corr(z_spec_corr>z_spec_max) = z_spec_max;
                img_corr(m,n,s,:) = z_spec_corr';
            end
        end
    end 
end

end