function indexes = StructArrayMatch(structArray,fieldName,desiredValue)

    fieldValues = getArrayFields(structArray,fieldName);
    indexes = find((fieldValues == desiredValue));

end
