function params = set_parameter(params, parameter_path, value)
%SET_PARAMETER Set a scalar value in a two-level parameter structure.

parts = strsplit(char(parameter_path), '.');
if numel(parts) ~= 2
    error('set_parameter:UnsupportedPath', ...
        'Expected a two-level path such as TSM.r_C, received %s.', parameter_path);
end
params.(parts{1}).(parts{2}) = value;
end

