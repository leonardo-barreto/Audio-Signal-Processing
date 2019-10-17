function extractedArray = getArrayFields(structArray,fieldName,arrayElements)

    extractedArray(1) = 0;

    if nargin > 2
        for index = 1:length(arrayElements)
            extractedArray(index) = getfield(structArray(arrayElements(index)),fieldName);
        end
    else
        extractedArray = arrayfun(@(singleStruct) getfield(singleStruct,fieldName),structArray);
    end

end