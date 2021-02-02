function [field,dims] = read_ft44_rfield(fid,ver,fieldname,dims,dims_unsure,extra_lines_after_label,soft_errors)
% field = read_ft44_rfield(fid,ver,fieldname,dims)
%
% Auxiliary routine to read real fields from fort.44 file
% 
%
% If you don't know a single dim, you can return it (no sanity check!)
% Ex = dims = [5,1,10], dims_unsure = [0,1,0]
% J.D. Lore, modified from W. Dekeyser
if nargin < 5
    dims_unsure = [];
end
if nargin < 6
    extra_lines_after_label = 0;
end
if nargin < 7
    soft_errors = 0;
end


if isempty(dims_unsure)
    FIND_DIM = 0;   
else
    FIND_DIM = 1;
    if sum(dims_unsure) ~= 1
        error('Do not understand, specify a single unknown dim')
    end    
    if dims(logical(dims_unsure)) ~= 1
        error('The given dim of the unknown dimension should be one')
    end
end
    
    
% Author: Wouter Dekeyser
% E-mail: wouter.dekeyser@kuleuven.be
% November 2016
% J.D. Lore modifications

% Verion 20160829: field label and size are specified in fort.44
% Do consistency check on the data
% if ver >= 20160829
    % Search the file until identifier 'fieldname' is found
    line = fgetl(fid);
    while isempty(strfind(line,fieldname))
        line = fgetl(fid);
        if line == -1
            if soft_errors
                fprintf('EOF reached without finding %s\n',fieldname)
                field = [];
                dims = [];
                return;
            else
                error(['EOF reached without finding ',fieldname,'.']);
            end
        end
    end
    for i = 1:extra_lines_after_label
        line = fgetl(fid);
        if line == -1
            if soft_errors
                fprintf('EOF reached without finding %s\n',fieldname)
                field = [];
                dims = [];                
                return;
            else
                error(['EOF reached without finding ',fieldname,'.']);
            end
        end        
    end
    
    
    % Consistency check: number of elements specified in the file should equal
    % prod(dims)
    % Handle cases with parentheses
    parPos = strfind(line,')');
    if ~isempty(parPos)
        numin = strread(line(parPos+1:end),'%*s %*s %d');
    else
        numin = strread(line,'%*s %*s %*s %*s %*s %*s %d');
    end
    if FIND_DIM
        dims(logical(dims_unsure)) = numin/prod(dims);
    end
    if numin ~= prod(dims)
        error('read_ft44_rfield: inconsistent number of input elements.');
    end
% end

% Read the data
field = fscanf(fid,'%e',prod(dims));
if (length(dims) > 1)
    field = reshape(field,dims);
end