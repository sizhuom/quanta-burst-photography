function f = defaultField(s, fn, dv)
%DEFAULTFIELD Get field value of a struct. If the field does not exist, a
%default value is returned

if isfield(s, fn)
    f = s.(fn);
else
    f = dv;
end

end

