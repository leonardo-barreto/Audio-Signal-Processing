function [path] = path_check(path)
%PATH_CHECK checks if a path exists and if not, it is created

if exist(path,'dir') ~= 7
    mkdir(path);
end

end