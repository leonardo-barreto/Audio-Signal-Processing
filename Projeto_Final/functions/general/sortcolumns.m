function matrixSorted = sortcolumns(matrix,row,direction)

switch nargin

    case 1
        matrixSorted = sortrows(transpose(matrix));
    case 2
        matrixSorted = sortrows(transpose(matrix),row);
    case 3
        if (~strcmp(direction,'ascend') && ~strcmp(direction,'descend'))
            error('Direction must be ''ascend'' or ''descend''');
        end
        matrixSorted = sortrows(transpose(matrix),row,direction);   
    otherwise
        error('Too many or too little arguments. 1 to 3 allowed.');

end

matrixSorted = matrixSorted';

end