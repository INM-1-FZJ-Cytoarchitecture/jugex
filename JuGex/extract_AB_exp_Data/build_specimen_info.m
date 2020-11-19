function [ specimen_info ] = build_specimen_info()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
str=urlread('http://api.brain-map.org/api/v2/data/query.json?criteria=model::Donor,rma::criteria,products[id$eq2],rma::include,age,rma::options[only$eq%27donors.id,donors.name,donors.race_only,donors.sex%27]');
json = parse_json(str);

for i=1:json.num_rows
    specimen_info(i).id=json.msg{1,i}.id;
    specimen_info(i).name=json.msg{1,i}.name;
    specimen_info(i).race=json.msg{1,i}.race_only;
    specimen_info(i).gender=json.msg{1,i}.sex;
    specimen_info(i).age=(json.msg{1,i}.age.days)/365;
end

end

