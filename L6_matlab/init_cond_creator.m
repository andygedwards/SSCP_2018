function init_cond_creator(filename, newfilename)
  data = load(filename);
  
  fid = fopen(newfilename, 'w');

style = '';
%for n = 1:(29 +6+ 4*2)
for n = 1:(29 +6+ 4*2 + 1) % to introduce new variable for IKCa
    style = [style ' %6.6e'];
end

%fprintf(fid, style, [data(end, 2:44)]);
%fprintf(fid, style, [data(end, 1:15) 0 data(end, 16:43)]);
fprintf(fid, style, [data(end, 2:45)]);

fclose(fid);
