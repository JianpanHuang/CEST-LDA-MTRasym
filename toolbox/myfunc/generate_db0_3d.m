function db0 = generate_db0_3d(img, offs, roi, method)
% FUNCTION:
%   To generate deltaB0 map using spline interpolation or Lorentzian fitting
%   Note: spline interpolation is faster but Lorentzian fitting is more accurate
% INPUT:
%   img: images used for deltaB0 map generation
%   offs: frequency offsets
%   roi: ROI for deltaB0 map generation
%   method: could be 'spline' or 'lorentz'
% OUTPUT:
%   db0: deltaB0 map
% REFERENCE
%   Kim M, et al., Magnetic Resonance in Medicine, 2009, 61(6): 1441-1450.
% AUTHOR:
%   Jianpan Huang, Email: jianpanhuang@outlook.com

%%
[xn, yn, sn, ~] = size(img);
if nargin < 3
    roi = ones(xn, yn, sn);
end
if nargin < 4
    method = 'spline';
end

% Generate deltaB0 map
db0 = zeros(xn, yn, sn);
offs_interp = min(offs(:)):0.01:max(offs(:));
h = waitbar(0,'Generating \DeltaB0 map, please wait >>>>>>');
for s = 1:sn
    if strcmp(method,'spline') == 1
        [row, col] = find(abs(roi(:,:,s))>0);
        for n = 1:length(row)
            z_spec = squeeze(img(row(n),col(n), s,:));
            z_spec_interp = spline(offs, z_spec, offs_interp);
            [~,ind] = min(z_spec_interp);
            shift_val = offs_interp(ind);
            if abs(shift_val) > 1 % if rNOE peak is lower than water peak, this checking procedure can fix it
                offs_interp_check = -abs(shift_val)+0.4:0.01:abs(shift_val)-0.4;
                z_spec_interp_check = spline(offs, z_spec, offs_interp_check);
                [~,ind_check] = min(z_spec_interp_check);
                ratio = ind_check/length(offs_interp_check);
                if ratio>0.2 && ratio<0.8
                    shift_val = offs_interp_check(ind_check);
                end
            end
            db0(row(n),col(n),s)= shift_val;
        end
    elseif strcmp(method,'lorentzian') == 1
        [row, col] = find(abs(roi(:,:,s))>0);
        h = waitbar(0,'Generating \DeltaB0 map, please wait >>>>>>');
        for n = 1:length(row)
            z_spec = squeeze (img(row(n),col(n), s,:));
            z_spec_interp = spline(offs,z_spec,offs_interp);
            [~,ind] = min(z_spec_interp);
            shift_val = offs_interp(ind);
            if abs(shift_val) > 1 % if rNOE peak is lower than water, this checking procedure can fix it
                offs_interp_check = -abs(shift_val)+0.4:0.01:abs(shift_val)-0.4;
                z_spec_interp_check = spline(offs, z_spec, offs_interp_check);
                [~,ind_check] = min(z_spec_interp_check);
                ratio = ind_check/length(offs_interp_check);
                if ratio>0.2 && ratio<0.8
                    shift_val = offs_interp_check(ind_check);
                end
            end
            % Lorentzian fit
            iv = [0.9,  0.1,   10,   shift_val];
            lb = [0.2,  0.01,  0.1,  shift_val-2];
            ub = [1,    1,     100,  shift_val+2];
            par = lsqcurvefit(@lorentz1pool, iv, offs, z_spec, lb, ub);
            db0(row(n),col(n), s) = par(4);
    %         z_spec_fit = lorentz1pool(par, offs);
    %         if mod(n, 100) == 0
    %             plot(offs, z_spec, 'o', offs, z_spec_fit, '-')
    %         end
        end
    end
    waitbar(s/sn,h);
end
delete(h);
end

function y = lorentz1pool(p, x)
% One-pool Lorentzian function
y = p(1) - p(2)*p(3)^2/4./(p(3)^2/4+(x-p(4)).^2) ;
end