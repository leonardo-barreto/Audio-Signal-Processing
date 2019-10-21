function [elements,indexes] = sortcolumns(matrix,row,direction)

switch nargin

    case 1
        [elements,indexes] = sortrows(transpose(matrix));
    case 2
        [elements,indexes] = sortrows(transpose(matrix),row);
    case 3
        if strcmp(direction,'ascend')
            [elements,indexes] = sortrows(transpose(matrix),row);
        elseif strcmp(direction,'descend')
            [elements,indexes] = sortrows(transpose(matrix),row);
            elements = flip(elements);
        else
            error('Direction must be ''ascend'' or ''descend''');
        end
    otherwise
        error('Too many or too little arguments. 1 to 3 allowed.');

end

elements = transpose(elements);
indexes = transpose(indexes);

end