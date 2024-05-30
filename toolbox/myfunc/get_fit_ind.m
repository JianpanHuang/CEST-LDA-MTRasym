function fit_ind = get_fit_ind(z_spec, offs, fit_rgc, fit_rge)
% FUNCTION:
%   To get index of Z-spectrum used of Lorentzian fitting
% INPUT:
%   z_spec: Z-spectrum
%   offs: frequency offsets
%   fit_rgc: center range, in ppm
%   fit_rge: edge range, in ppm
% OUTPUT:
%   fit_ind: fitting images
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
offs_ind = 1:length(offs);
% Edge, positive side
offs_flag = offs>fit_rge(2);
if sum(offs_flag) ~= 0
    fit_ind_pos = offs_ind(offs_flag);
    z_spec_pos = z_spec(offs_flag);
    mean_temp = mean(z_spec_pos);
    std_temp = std(z_spec_pos);
    fit_ind_pos(z_spec_pos>mean_temp+std_temp | z_spec_pos<mean_temp-std_temp ) = []; 
else
    fit_ind_pos = [];
end
% Edge, negative side
offs_flag = offs<fit_rge(1);
if sum(offs_flag) ~= 0
    fit_ind_neg = offs_ind(offs_flag);
    z_spec_neg = z_spec(offs_flag);
    mean_temp = mean(abs(z_spec_neg));
    std_temp = std(abs(z_spec_neg));
    fit_ind_neg(abs(z_spec_neg)>mean_temp+std_temp | abs(z_spec_neg)<mean_temp-std_temp ) = []; 
else
    fit_ind_neg = [];
end
% Center
fit_ind_cent = offs_ind(offs>fit_rgc(1) & offs<fit_rgc(2));

% Check existence
if isempty(fit_ind_pos) && isempty(fit_ind_neg)
    fit_ind_pos = offs_ind(end); % make sure at least one point at edges is included
end
if isempty(fit_ind_cent)
    [~, min_ind]=min(z_spec);
    fit_ind_cent = offs_ind(min_ind-2:min_ind+2); % make sure at least five point at center is included
end

% Combine
fit_ind = [fit_ind_neg, fit_ind_cent, fit_ind_pos];

end