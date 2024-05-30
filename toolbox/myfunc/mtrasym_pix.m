function [img_mtrasym, img_mtrasymrex] = mtrasym_pix(img, offs, roi)
% FUNCTION:
%   To do MTRasym and MTRasymRex pixel by pixel
%   Note: both conventional and inverse Z-spectrum analysis methods are included
% INPUT:
%   img: Z-spectrum images
%   offs: frequency offsets
%   roi: ROI
% OUTPUT:
%   img_mtrasym: MTRasym images
%   img_mtrasymrex: MTRasymRex images
% REFERENCE
%   MTRasym: Zhou J, et al. Prog Nucl Magn Reson Spectrosc 2006;48(23):109-136.
%   Inverse Z-spectrum analysis: Zaiss M, et al., NMR in biomedicine, 2014, 27(3): 240-252.
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
if nargin < 4
    roi = ones(size(img));
end

% MTRasym and MTRasymRex
sz = size(img);
offs_pos = offs(offs>=0);
offs_max = max(abs(offs));
offs_pos_interp = 0:0.01:offs_max;
offs_neg_interp = 0:-0.01:-offs_max;
img_mtrasym = zeros(sz(1),sz(2),length(offs_pos));
img_mtrasymrex = zeros(sz(1),sz(2),length(offs_pos));
for m = 1:sz(1)
    for n = 1:sz(2)
        if roi(m,n) == 1
            z_spec_temp = squeeze(img(m,n,:));
            z_spec_pos_interp  = interp1(offs,z_spec_temp,offs_pos_interp,'spline');
            z_spec_neg_interp  = interp1(offs,z_spec_temp,offs_neg_interp,'spline');
            % MTRasym
            mtrasym_temp = z_spec_neg_interp-z_spec_pos_interp;
            mtrasym = interp1(offs_pos_interp,mtrasym_temp,offs_pos,'spline');
            img_mtrasym(m,n,:) = mtrasym;
            % MTRasymRex
            mtrasym_rex_temp = 1./(z_spec_pos_interp+eps)-1./(z_spec_neg_interp+eps);
            mtrasym_rex = interp1(offs_pos_interp,mtrasym_rex_temp,offs_pos,'spline');
            img_mtrasymrex(m,n,:) = mtrasym_rex;
        end
    end
end