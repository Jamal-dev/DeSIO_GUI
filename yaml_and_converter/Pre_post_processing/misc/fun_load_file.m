function out = fun_load_file(fullFileName)

if exist(fullFileName, 'file')
    out = load(fullFileName);
else
  out = [];
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName);
end

return 
