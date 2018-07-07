function[] = struct2vars(varstruct)
if isstruct(varstruct)
    varnames = fieldnames(varstruct);
    
    for ii = 1:length(varnames)
        assignin('caller',varnames{ii},varstruct.(varnames{ii}))
    end
else
    error('Input is not a struct.')
end
end
