function [structto] = copy_class_fields_to_struct(classfrom, structto, prefix)
%COPY_CLASS_FIELDS_TO_STRUCT Summary of this function goes here
%   Detailed explanation goes here
arguments
    classfrom 
    structto

    % Optional
    prefix = ""
end

if ismember("fields_to_copy", fieldnames(classfrom))
    for f = 1:numel(classfrom.fields_to_copy)
        field_name = classfrom.fields_to_copy(f);
        structto.(prefix + field_name) = classfrom.(field_name);
    end
else
    for n = 1:numel(fieldnames(classfrom))
        fn = fieldnames(classfrom);
        if fn(n) ~= "Name"
            structto.(prefix + fn{n}) = classfrom.(fn{n});
        end
    end
end

end

