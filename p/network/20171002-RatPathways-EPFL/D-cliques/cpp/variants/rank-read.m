
%%

input_file = './test_J2I_example.txt';
fid = fopen(input_file);
I = [];
J = [];
j = 0;
while 1
    line = fgetl(fid);
    if (~ischar(line)); break; end
    j = j + 1;
    
    I = [I, 1 + str2double(strsplit(line))];
    J = [J, j * ones(1, length(I) - length(J))];
    
    if (mod(j, 100) == 0); disp(j); end
end
fclose(fid);

M = sparse(I, J, ones(size(I)));

spy(M)

%%


