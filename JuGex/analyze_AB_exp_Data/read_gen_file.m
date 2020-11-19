function [probe_id,gene_symbol,entrez_id] = read_gen_file(file2read)
%read gen file
filename = file2read;
delimiter = ',';
startRow = 2;
formatSpec = '%f%s%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
probe_id = dataArray{:, 1};
gene_symbol = dataArray{:, 2};
entrez_id = dataArray{:, 3};
clearvars filename delimiter startRow formatSpec fileID dataArray ans;
end