function [flag, idx] = my_isfield(structure, field_name) 
% FUNCTION:
% 	To search if 'field_name' exist in 'structure'
% INPUT:
%   structure: structure to search from
%   field_name: name of field to be searched
% OUTPUT:
%   flag: 1 for yes, 0 for no
% AUTHOR: 
%   Jianpan huang, Email: jianpanhuang@outlook.com

%%
flag = 0; 
all_field = fieldnames(structure(1)); 
for i=1:length(all_field) 
    if(strcmp(all_field{i}, strtrim(field_name))) 
        flag = 1; 
        idx = i;
        return; 
    end
end
end