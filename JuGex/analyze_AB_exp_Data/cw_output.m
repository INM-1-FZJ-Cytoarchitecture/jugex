function cw_output(ind_sig_gens,probe_id,entrez_id,gene_symbol,p,str_title,table_header,corr_type,pused,savestr,search_mode,n_rep)

if isempty(ind_sig_gens)
    disp([char(10) char(10) str_title]);
    disp('no significant results');
else
    probs=probe_id(ind_sig_gens);
    entrezs=entrez_id(ind_sig_gens);
    symbols=gene_symbol(ind_sig_gens);
    adj_ps=p(ind_sig_gens);
    
    disp([char(10) char(10) str_title]);
    format long;
    if (size(ind_sig_gens,1)<size(ind_sig_gens,2))
        tmp_ind_sig_gens=ind_sig_gens';
    else
        tmp_ind_sig_gens=ind_sig_gens;
    end
    if (size(adj_ps,1)<size(adj_ps,2))
        tmp_adj_ps=adj_ps';
    else
        tmp_adj_ps=adj_ps;
    end
    if (size(probs,1)<size(probs,2))
        tmp_probs=probs';
    else
        tmp_probs=probs;
    end
    if (size(entrezs,1)<size(entrezs,2))
        tmp_entrezs=entrezs';
    else
        tmp_entrezs=entrezs;
    end
    if (size(symbols,1)<size(symbols,2))
        tmp_symbols=symbols';
    else
        tmp_symbols=symbols;
    end
    
    
    if search_mode==1 %genes
        try
            if isdir(['.' filesep 'output' filesep 'analyzed_data' filesep 'all-probes mode'])~=1
                mkdir(['.' filesep 'output' filesep 'analyzed_data' filesep 'all-probes mode']);
            end
        end
        %   table_header={'Index' 'Entrez_ID' 'Gene_Symbol' 'p_FWE_corrected'};
        T = table(tmp_ind_sig_gens,tmp_entrezs,tmp_symbols, tmp_adj_ps,'VariableNames',table_header')
        writetable(T,['.' filesep 'output' filesep 'analyzed_data' filesep 'all-probes mode' filesep savestr '_' corr_type '_' pused '_n_rep_' num2str(n_rep) '.csv'])
    elseif search_mode==2 %probes
        try
            if isdir(['.' filesep 'output' filesep 'analyzed_data' filesep 'single-probe mode'])~=1
                mkdir(['.' filesep 'output' filesep 'analyzed_data' filesep 'single-probe mode']);
            end
        end
        if size(tmp_probs,1)>1
            for i=1:size(tmp_probs,1)  % tmp_probs(2,1)
                url=['http://api.brain-map.org/api/v2/data/Probe/' num2str(tmp_probs(i,1)) '.json?include=predicted_sequence'];
                status=0;
                while status<1
                    try
                        [str,status] = urlread(url);
                    catch
                        disp('lost connection. I start a new attempt!');
                    end
                end
                json = parse_json(str);
                probe_name_cell{i}=json.msg{1,1}.name;
            end
            probe_name_cell=probe_name_cell';
        else
            url=['http://api.brain-map.org/api/v2/data/Probe/' num2str(tmp_probs) '.json?include=predicted_sequence'];
            status=0;
            while status<1
                try
                    [str,status] = urlread(url);
                catch
                    disp('lost connection. I start a new attempt!');
                end
            end
            json = parse_json(str);
            probe_name_cell{1}=json.msg{1,1}.name;
        end
        % probe label =     json.msg{1,1}.name
        %   table_header={'Index' 'Probe_name' 'Entrez_ID' 'Gene_Symbol' 'p_FWE_corrected'};
        T = table(tmp_ind_sig_gens,probe_name_cell,tmp_entrezs,tmp_symbols,tmp_adj_ps,'VariableNames',table_header')
        writetable(T,['.' filesep 'output' filesep 'analyzed_data' filesep 'single-probe mode' filesep savestr '_' corr_type '_' pused '_n_rep_' num2str(n_rep) '.csv'])
    end
    %writetable(T,['.' filesep 'output' filesep 'analyzed_data' filesep savestr '_' corr_type '_' pused '.xls'])
end
