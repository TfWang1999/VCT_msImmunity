function value = get_parameter(params, parameter_path)
%GET_PARAMETER Read a scalar value from a nested parameter structure.

parts = strsplit(char(parameter_path), '.');
value = params;
for i = 1:numel(parts)
    value = value.(parts{i});
end
end

