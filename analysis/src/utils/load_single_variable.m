function value = load_single_variable(file_path)
%LOAD_SINGLE_VARIABLE Load the only top-level variable stored in a MAT file.

contents = load(file_path);
names = fieldnames(contents);
if numel(names) ~= 1
    error('load_single_variable:VariableCount', ...
        'Expected one variable in %s, found %d.', file_path, numel(names));
end
value = contents.(names{1});
end

