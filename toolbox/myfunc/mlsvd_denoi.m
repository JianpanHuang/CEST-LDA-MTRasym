function img_denoi = mlsvd_denoi(img, trunc_ind)
% FUNCTION
%   To denoise input images
% INPUT
%   img: images
%   trunc_ind: truncation index
% OUTPUT
%   img_denoi: denoising images
% AUTHOR
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%% 
% Define truncation index if it is not given
if nargin < 2
    % 1st method: roughly determine truncation index based image size
    [xn, yn, wn] = size(img);
    trunc_ind = round([xn/3, yn/3, wn/4]);
%     % 2nd method: determine truncation index based on Malinowskis, Nelson and Median criteria
%     % See reference: Huang J. et al., NMR in Biomedicine, 2022, 35(3): e4640.
%     [~, ~, svall] = mlsvd(img);
%     sv{1} = svall{1}/max(svall{1});
%     sv{2} = svall{2}/max(svall{2}); 
%     sv{3} = svall{3}/max(svall{3}); 
%     [mal_ind(1), nel_ind(1), med_ind(1)] = determ_trunc_index(sv{1});
%     [mal_ind(2), nel_ind(2), med_ind(2)] = determ_trunc_index(sv{2});
%     [mal_ind(3), nel_ind(3), med_ind(3)] = determ_trunc_index(sv{3});
%     trunc_ind = [med_ind(1), med_ind(2), nel_ind(3)];
end

% Denoise images
[ue, se, ~] = mlsvd(img, trunc_ind);
img_denoi = lmlragen(ue, se);
end