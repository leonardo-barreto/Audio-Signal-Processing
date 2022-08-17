function [elements,indexes] = sortcolumns2018b(matrix,row,direction)

switch nargin

    case 1
        [elements,indexes] = sortrows(transpose(matrix));
    case 2
        [elements,indexes] = sortrows(transpose(matrix),row);
    case 3
        [elements,indexes] = sortrows(transpose(matrix),row,direction);

end

elements = transpose(elements);
indexes = transpose(indexes);

end