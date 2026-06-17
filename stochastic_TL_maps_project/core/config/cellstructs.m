function s = cellstructs(c)
% Convert a cell array of structs (possibly with different fields) to a struct array.
% Missing fields are filled with [].
    all_fields = {};
    for i = 1:numel(c)
        all_fields = union(all_fields, fieldnames(c{i}));
    end
    s = struct();
    for i = 1:numel(c)
        for f = all_fields'
            if isfield(c{i}, f{1})
                s(i).(f{1}) = c{i}.(f{1});
            else
                s(i).(f{1}) = [];
            end
        end
    end
end
