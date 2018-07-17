function[] = struct2vars(vars,varstruct)
if isstruct(varstruct)
    varnames = fieldnames(varstruct);
    finalvars = {varnames{ismember(varnames,vars)}}; 
    for ii = 1:length(finalvars)
        assignin('caller',finalvars{ii},varstruct.(finalvars{ii}))
    end
else
    error('Input is not a struct.')
end
end
