function db0 = generate_load_db0_3d(path, db0_name, img, offs, roi, method)
% FUNCTION:
%   To generate or load deltaB0 map
% INPUT:
%   path: file path
%   db0_name: name of deltaB0 map
%   img: images used for deltaB0 map generation
%   offs: frequency offsets
%   roi: ROI for deltaB0 map generation
%   method: could be 'spline' or 'lorentzian'
% OUTPUT:
%   db0: deltaB0 map
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
if ~exist([path, filesep, db0_name, '.mat'], 'file') == 1
    db0 =  generate_db0_3d(img, offs, roi, method);
    save([path, filesep, db0_name, '.mat'], 'db0') 
else
    load([path, filesep, db0_name, '.mat'], 'db0');
end
end