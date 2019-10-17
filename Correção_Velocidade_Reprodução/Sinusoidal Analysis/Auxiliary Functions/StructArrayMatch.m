function indexes = structArrayMatch(structArray,fieldName,desiredValue)

    fieldValues = getArrayFields(structArray,fieldName);
    indexes = find((fieldValues == desiredValue));

end
