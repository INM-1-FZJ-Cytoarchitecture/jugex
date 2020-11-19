function [ combined_zscores,area1_zscores,area2_zscores,unique_entrez_id  ] = switch2gensymbol_lvl( entrez_id,combined_zscores,size_area1_zscores,size_area2_zscores )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

unique_entrez_id=unique(entrez_id);
winsorzed_mean_zscores=zeros(size(combined_zscores,1),size(unique_entrez_id,1));

for i=1:size(unique_entrez_id,1)
    tmp=combined_zscores(:,find(entrez_id==unique_entrez_id(i))); %index aller values einer unique entrez_id. davon dann die zscores
    for j=1:size(combined_zscores,1)
        %mean(winsor(tmp(j,:),[10 90]));
        if size(tmp,2)==1
            %disp(['only one probe ' unique_entrez_id(i)]);
            winsorzed_mean_zscores(j,i)=tmp(j,1);
        else
        winsorzed_mean_zscores(j,i)=mean(winsor(tmp(j,:),[10 90]));
        end
    end
end

% zurück schreiben von Area1 & Area2 zscores um nur noch winsorzed means zu haben
combined_zscores=winsorzed_mean_zscores;
area1_zscores=winsorzed_mean_zscores(1:size_area1_zscores,:);
area2_zscores=winsorzed_mean_zscores(size_area1_zscores+1:size_area1_zscores+size_area2_zscores,:);


end

