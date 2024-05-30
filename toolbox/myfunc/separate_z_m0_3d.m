function [img, offs, img_m0] = separate_z_m0_3d(img_all, offs_all, offs_m0)
% FUNCTION:
%	To separate M0 image from Z-spectrum images
% INPUT:
%	img_all: all images
%	offs_all: all offsets
%	offs_m0: M0 offset
% OUTPUT:
% 	img: Z-spectrum images
%	offs: Z-spectrum offsets
%   img_m0: M0 image
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
% Set M0 offset as 99 ppm if it is not given, 
if nargin < 3
    offs_m0 = 99;
end

% Check size of images and offsets
img_all_dim = ndims(img_all);
if size(img_all, img_all_dim) ~= length(offs_all)
    errordlg('Image size and offset size do not match!','Data error');
    error('Image size and offset size do not match!');
end

% Separate M0 image from Z-spectrum images,
% and rearrange offsets and Z-spectrum images in ascending order
index = abs(offs_all) > abs(offs_m0)-1;
if sum(index) ~= 0
    img_m0_all = img_all(:,:,:,index);
    m0_dim = ndims(img_m0_all);
    if m0_dim == 4
        img_m0 = img_m0_all(:,:,:,end);
    else
        img_m0 = img_m0_all;
    end
    offs_temp = offs_all;
    offs_temp(index) = [];
    img_temp = img_all;
    img_temp(:,:,:,index) = [];
    [offs,pos] = unique(offs_temp);
    img = img_temp(:,:,:,pos);
else
    [offs, pos] = unique(offs_all);
    img = img_all(:,:,:,pos);
    % Take the image with the largest offset as M0 image if there is no M0 image
    img_m0 = img_all(:,:,:,end);
end
end