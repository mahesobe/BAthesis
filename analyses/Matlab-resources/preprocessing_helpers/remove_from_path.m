function [] = fv_remove_from_path(target)
tmpPath = path;
strIdx = findstr(path,target);

for k = strIdx(length(strIdx):-1:1)
   firstIdx = findstr(tmpPath(1:k),':');
   lastIdx = findstr(tmpPath(k:end),':');
   if isempty(firstIdx)||isempty(lastIdx)
       continue;
   end
   tmpPath(firstIdx(end)+1:k+lastIdx(1)-1) = [];    
end

path(tmpPath)
