function out=parse_aba_probes_file(filename,gene_liste_output_name)
% parse_aba_probes_file(filename,gene_liste_output_name)
%
%   read gene symbols list (text file) and extracts Allen Brain probe ids
%   It is checked whether genes were analyzed by Allen Brain, if not the
%   missing gene symbols are listed in the console
%
%   filename:  Input text file, one gene symbol per row
%   gene_liste_output_name: name of formated gene list. This file should be
%   copied to the gene list folder
% 
% 
% 
% Input file structur:
% DRD1
% CHRM1
% CHRM2
% CHRM3
% CHRNA1
% 
% Output file structure:
% probe_id,gene_symbol,entrez_id
% 1022773,"DRD1",1812
% 1022768,"DRD1",1812
% 1022767,"DRD1",1812
% 1022757,"DRD1",1812
% 1022756,"DRD1",1812
% 1022752,"DRD1",1812
% 1022751,"DRD1",1812
% 1022747,"DRD1",1812
% 1022746,"DRD1",1812
% 1022745,"DRD1",1812
% 1022744,"DRD1",1812
% 1022743,"DRD1",1812
% 1022742,"DRD1",1812
% 1022741,"DRD1",1812
% 1022740,"DRD1",1812
% 1022739,"DRD1",1812
% 1022738,"DRD1",1812
% 1022737,"DRD1",1812
% 1022736,"DRD1",1812
% 1022735,"DRD1",1812
% 1058251,"CHRM1",1128
% 1058242,"CHRM1",1128
% ....
%   
%   example:
%   parse_aba_probes_file('text_file.txt','gene_liste_output_name.csv');
%   
%   text_file.txt STRUCTURE
%   no header line
%   gene_symbols as string
%   only line break after each row, no comma
%
%
%   read complete file into table and saves gene list as csv file


formatspec=('%d  %q  %d  %q   %q  %d  %q');
T = readtable('Probes.csv','Delimiter',',', 'Format',formatspec);
search_symbols = read_txt_gene_file(filename);
search_symbols=search_symbols';
%% search_symbols cell array MDD associated genes
%search_symbols={'AAK1','AKT3','CAMK1G','CAMK4','CDH4','DLG2','EIF4B','GALNT13','LRRFIP1','MPPED2','MTG2','MTHFD1L','NTRK3','RIMBP2','RYR2','TRAF3','VEGFA','WWOX'};
%% rows = n x 1 logical array => 1 at positions where gene_symbol was found
result_t=table;
for i=1:size(search_symbols,2)
    rows = strcmp(T.gene_symbol,search_symbols{i});
    %% output vars of T
    vars = {'probe_id','gene_symbol','entrez_id'};
    %% values of vars of T where gene_symbol was found 
    result_t = [result_t;T(rows,vars)];
    
end

if size(unique(result_t.gene_symbol),1)~=size(search_symbols,2)
disp('Following gene symbols were not found:');
c = setdiff(search_symbols,unique(result_t.gene_symbol)')

end
% save as csv 
writetable(result_t,gene_liste_output_name,'Delimiter',',','QuoteStrings',true);
out='done';
