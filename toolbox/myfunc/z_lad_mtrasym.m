function z_res = z_lad_mtrasym(path, img, offs, par, roi)
% FUNCTION:
%   To analyze Z-spectra using Lorentzian difference analysis (LDA) and MTRasym
% INPUT:
%   path: file path
%   img: Z-spectrum images
%   offs: frequency offsets
%   par: other parameters
%   roi: ROI
% OUTPUT:
%   cest_res: all cest analysis results
% REFERENCE:
%   MTRasym: Zhou J, et al. Prog Nucl Magn Reson Spectrosc 2006;48(23):109-136.
%   LDA: Zaiß M, et al., Journal of magnetic resonance, 2011, 211(2): 149-155.
%        Jones C K, et al., Magnetic resonance in medicine, 2012, 67(6): 1579-1589.
%   Inverse Z-spectrum analysis: Zaiss M, et al., NMR in biomedicine, 2014, 27(3): 240-252.
%   AREX method: Zaiss M, et al., Neuroimage, 2015, 112: 180-188.
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
% Do Lorentzian fitting or load fitted Lorentzian data
if nargin < 3
    roi = ones(xn, yn);
end
if my_isfield(par,'t1map')
    t1map = par.t1map;
    t1map(roi==1&t1map==0) = mean(t1map(t1map>0));
    r1map = 1./(t1map+eps);
else
    r1map = ones(size(img,1),size(img,2)); % If no t1 map, all t1 are set to 1 s.
end

if my_isfield(par,'file_name')
    save_path = [path,filesep,par.file_name,'_z_res.mat'];
else
    save_path = [path,filesep,'z_res.mat'];
end

% Analysis
list = {'Import the fitting data','Do Lorentzian fitting'};
[list_num,~] = listdlg('ListString',list);
if list_num == 2
    [img_lf, img_lda, img_ldarex] = lorentz_fit_pix(img, offs, par, roi);
    img_ldaarex = img_ldarex.*repmat(r1map,[1,1,size(img_ldarex,3)]);
    [img_mtrasym, img_mtrasymrex] = mtrasym_pix(img, offs, roi);
    img_mtrasymarex = img_mtrasymrex.*repmat(r1map,[1,1,size(img_mtrasymrex,3)]);
    save(save_path,'img_lf','img_lda','img_ldarex','img_ldaarex',...
                                       'img_mtrasym','img_mtrasymrex','img_mtrasymarex');
else
    if ~exist(save_path,'file') ==1
        [img_lf, img_lda, img_ldarex] = lorentz_fit_pix(img, offs, par, roi);
        img_ldaarex = img_ldarex.*r1map;
        [img_mtrasym, img_mtrasymrex] = mtrasym_pix(img, offs, roi);
        img_mtrasymarex = img_mtrasymrex.*r1map;
        save(save_path,'img_lf','img_lda','img_ldarex','img_ldaarex',...
                                           'img_mtrasym','img_mtrasymrex','img_mtrasymarex');
    else
        load(save_path,'img_lf','img_lda','img_ldarex','img_ldaarex',...
                                           'img_mtrasym','img_mtrasymrex','img_mtrasymarex');
    end
end
z_res.img_lf = img_lf;
z_res.img_lda = img_lda;
z_res.img_ldarex = img_ldarex;
z_res.img_ldaarex = img_ldaarex;
z_res.img_mtrasym = img_mtrasym;
z_res.img_mtrasymrex = img_mtrasymrex;
z_res.img_mtrasymarex = img_mtrasymarex;
end