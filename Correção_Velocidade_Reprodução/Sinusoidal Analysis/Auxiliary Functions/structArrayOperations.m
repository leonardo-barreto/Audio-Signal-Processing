function indexes = structArrayOperations(structArray,fieldName,operation,desiredValue)

    fieldValues = getArrayFields(structArray,fieldName);

    if (~strcmp(fieldName,'status'))
        switch operation
            case '=='
                indexes = find(fieldValues == desiredValue);
            case '>'
                indexes = find(fieldValues > desiredValue);
            case '>='
                indexes = find(fieldValues >= desiredValue);
            case '<'
                indexes = find(fieldValues < desiredValue);
            case '<='
                indexes = find(fieldValues <= desiredValue);
            case '~='
                indexes = find(fieldValues ~= desiredValue);
            otherwise
                error('Insert a valid operation (as a character): ==,>,<,~=');
        end
    else
        switch desiredValue
            case 'inactive'
                indexes = find(fieldValues == 0);
            case 'active'
                indexes = find(fieldValues == 1);
            case 'asleep'
                indexes = find(fieldValues == 2);
        end
    end
    
end
