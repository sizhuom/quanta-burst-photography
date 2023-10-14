function success = checkStructFields(s, ref, fields)
%CHECKSTRUCTFIELDS Check if the fields of s are consistent with ref

if nargin < 3
    fields = fieldnames(ref);
end
for i = 1:numel(fields)
    if ~(isfield(s, fields{i}) && isequal(s.(fields{i}), ref.(fields{i})))
        success = false;
        return
    end
end
success = true;
end

