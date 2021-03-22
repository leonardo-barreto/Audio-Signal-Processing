function [path] = path_check_del(path)
%PATH_CHECK checks if a path exists. If yes, the directory is erased. If no, it is created.

if exist(path,'dir') == 7
    rmdir(path,'s');
end

mkdir(path);

end