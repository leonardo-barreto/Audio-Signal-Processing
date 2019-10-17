function indexes = structArrayOperations(structArray,fieldName,operation,desiredValue)

    fieldValues = getArrayFields(structArray,fieldName);

    switch operation
        case '=='
            indexes = find((fieldValues == desiredValue));
        case '>'
            indexes = find((fieldValues > desiredValue));
        case '<'
            indexes = find((fieldValues < desiredValue));
        case '~='
            indexes = find((fieldValues ~= desiredValue));
        otherwise
            error('Insert a valid operation (as a character): ==,>,<,~=');
    end
end
