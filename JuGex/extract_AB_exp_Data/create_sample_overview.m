function [ overview ] = create_sample_overview( input_args )
%create_sample_overview Summary of this function goes here
%   Detailed explanation goes here
areas=unique({input_args.name});

for i=1:size(input_args,2)
    area=input_args(i).name;
    IndexC = strfind(areas,area);
    Index = find(not(cellfun('isempty', IndexC)));
    
    overview(Index).name=area;
    rep_specimen = strrep(input_args(i).specimen, '.', '_');
    overview(Index).(rep_specimen)=size(input_args(i).validated_zscores,1);
    
    if isfield(overview(Index), 'n_samples')
        if isempty(overview(Index).n_samples)
            overview(Index).n_samples=0;
        end
        overview(Index).n_samples= overview(Index).n_samples+size(input_args(i).validated_zscores,1);
    else
        overview(Index).n_samples= size(input_args(i).validated_zscores,1);
    end
    
end
overview = orderfields(overview);
rep=fliplr([1:size(fieldnames(overview),1)]);
overview =orderfields(overview,rep);
end

