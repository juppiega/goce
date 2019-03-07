function n = linecount(fid)
frewind(fid);
n = 0;
tline = fgetl(fid);
while ischar(tline)
  tline = fgetl(fid);
  n = n+1;
end
frewind(fid);

end