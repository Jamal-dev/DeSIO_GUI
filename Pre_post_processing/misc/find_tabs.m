function [ntabs] = find_tabs(line)
%
% This function detects any "tab" characters in a character string and
% returns logical TRUE if tabs are present, or logical FALSE if tabs are
% not present. Optionally returns the positions of the tab(s) in the line.
%
% USAGE
% [x,tabpositions]=istab(line)
%
% INPUT
% line - [1 x n] character string with or without tabs
%
% OUTPUT
% x - logical TRUE or FALSE indicating presence of tab(s) in line.
% tabpositions - [1 x ntabs] vector showing positions of tabs in line.
%
% line must be a character string or empty
if ~ischar(line) & ~isempty(line)
    disp(' ')
    disp(' Input must be a character string.')
    return
end
% line must contain only one line
if size(line,1)>1
    disp(' ')
    error(' Character string must be no more than one line.' )
    return
end
% Put logical 1 where tab exists, logical 0 where it does not in string.
% Could use x=(line==9), but the following eliminates the case where this
% designator for the tab character is changed in the future.
x=(line==sprintf('\t'));
% Find the positions in line where tabs exist.
tabpositions=find(x);
if isempty(tabpositions), tabpositions=[]; end
% Return only logical 1 or 0, however many tabs there are in line.
x=x(find(x));
if ~isempty(x), x=x(1); end
if isempty(x), x=logical(0); end
% Display the logical result for tab existence.
if x 1==1; end
if ~x 1==2; end
ntabs = length(tabpositions);
return