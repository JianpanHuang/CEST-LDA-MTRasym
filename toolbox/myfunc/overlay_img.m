function overlay_img(img,img_bcg,roi,cax_rg,cmap,horiz_pos,fontsz,fontname)
% FUNCTION:
%   To overlay image on top of background image
% INPUT:
%	img: image on top
%   img_bcg: background image
%   roi: ROI mask
%   cax_rg: colorbar range
%   cmap: colormap type, e.g. parula, jet, hot, etc
%   horiz_pos: horizontal position of image, adjust this value if there is 
%              any horizontal shift because of different screen size
%   fontsz: font size of colorbar
% AUTHOR
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
if nargin < 3
    roi = zeros(size(img));
    roi(abs(roi(img))>0)=1;
end
if nargin < 4
    cax_rg = [min(img(:)),max(img(:))];
end
if nargin < 5
    cmap = jet;
end
if nargin < 6
    horiz_pos = 0.062;
end
if nargin < 7
    fontsz = 16;
end
if nargin < 8
    fontname = 'Arial';
end
if size(cax_rg) == 1
    if cax_rg < max(img(:))
        cax_rg(2) = max(img(:));
    else
        cax_rg(2) = cax_rg(1);
        cax_rg(1) = min(img(:));
    end
end
roi = uint8(roi(:,:,[1,1,1]));
img_bcg = imresize(img_bcg,size(img));

% Show image
scr_sz = get(0,'ScreenSize');
axes('position',[0.1,0.1,0.8,0.8]);
set(gcf,'Position',[scr_sz(3)*0.3, scr_sz(4)*0.3, scr_sz(3)*0.4, scr_sz(4)*0.5], 'color','w');
imshow(img,[], 'InitialMagnification','fit');colormap(cmap);
colorbar('FontSize',fontsz, 'FontName', fontname, 'FontWeight','bold' ,'linewidth',2);
caxis(cax_rg);
% Scale image
img = img - cax_rg(1);
img = img/(cax_rg(2)-cax_rg(1)); % Scale image to colorbar scale
len_cmap = length(cmap);
img_rgb = ind2rgb(gray2ind(img,len_cmap),cmap);
img_rgb = uint8(img_rgb*255);
% Scale background image
img_bcg = img_bcg - min(img_bcg(:));
img_bcg = img_bcg/max(img_bcg(:));
img_bcg_rgb = uint8(img_bcg(:,:,[1,1,1])*255);
% Overlay image
img_ol = img_rgb.*roi+img_bcg_rgb.*(1-roi);
% Show overlayed image
axes('position',[horiz_pos,0.1,0.8,0.8]);
imshow(img_ol,[], 'InitialMagnification','fit');
set(gcf,'Position',[scr_sz(3)*0.3, scr_sz(4)*0.3, scr_sz(3)*0.4, scr_sz(4)*0.5], 'color','w');
colorbar off;
end