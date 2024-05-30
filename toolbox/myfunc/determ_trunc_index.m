function [mal_ind, nel_ind, med_ind] = determ_trunc_index(sv)
% FUNCTION
%   To adaptively determine truncation index for MLSVD denoising based on Malinowskis, Nelson and Median criteria
% INPUT:
%   sv: Singular values
% OUTPUT:
%	malInd: The truncation index suggested by Malinowskis criterion
%	nelInd: The truncation index suggested by Nelson criterion
%	medInd: The truncation index suggested by Median criterion
% REFERENCES:
%   Breitling J. et al., NMR in Biomedicine, 2019, 32(11): e4133.
%   Huang J. et al., NMR in Biomedicine, 2022, 35(3): e4640.
% AUTHORS: 
%   Jianpan Huang, Email: jianpanhuang@outlook.com
%   Original version written by Johannes Breitling is available at https://github.com/jbreitling/CEST-AdaptiveDenoising

%%
d = sort(sv,'descend');
dm = diag(d);

% Malinowskis criterion. Malinowski, E. R. et al., Anal Chem 1977;49:612-617.
k_ind = zeros(length(d)-1,1);  
for mm = 1: length(d)-1
    k_ind(mm) = (sum(diag(dm(mm:end, mm:end))./(length(d)-mm)))^0.5 / (length(d)-mm)^2;
end
[~, mal_ind] = min(k_ind);

% Nelson criterion. Nelson, L.R. et al., J Educ Res Meas 2005;3:1-17.
for nn = 1:length(d)-1
    xx = nn:length(d);
    d_sub = d(nn:end);
    p = polyfit(xx', d_sub, 1);
    d_fit = p(1)*xx' + p(2);
    rsq(nn) = 1 - sum((d_sub - d_fit).^2)/sum((d_sub - mean(d_sub)).^2);
end
nel_ind = find(rsq > 0.80, 1, 'first'); % Find index with R^2 > 0.80.

% Median criteria. Manjón, J.V. et al., Med Image Anal. 2015;22:35-47.
beta = 1.29;
d_sub = d(sqrt(d) < 2*median(sqrt(d)));
med_ind = find(sqrt(d) >= beta*sqrt(median(d_sub)),1, 'last');
end