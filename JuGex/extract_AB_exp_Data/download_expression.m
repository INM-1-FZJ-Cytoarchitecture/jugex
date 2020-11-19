% Download expression data for a given set of probe_ids.
function [samples,explevels,zscores,probes] = download_expression(probe_ids, donor_ids,search_mode,struct_id,disp_url)

% Copyright 2013 Allen Institute for Brain Science
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% Convert the probe ids into a comma-separated string.
nprobes = numel(probe_ids);
if (nprobes > 1)
    probe_str = num2str(probe_ids','%d,');
    probe_str = probe_str(1:end-1);   
else
    probe_str = num2str(probe_ids','%d');
end

% Convert the donor ids into a comma-separated string.
ndonors = numel(donor_ids);
if (ndonors > 1)
    donor_str = num2str(donor_ids','%d,');
    donor_str = donor_str(1:end-1);
else
    donor_str = num2str(donor_ids','%d');
end
     
if search_mode==2
    structure_str = num2str(struct_id);
    structures = sprintf('[structures$eq%s]',structure_str);
end
    


% Make the request and parse the results as JSON.
service = get_api_path();
request = '/query.json?criteria=service::human_microarray_expression';
probes = sprintf('[probes$in%s]',probe_str);
donors = sprintf('[donors$eq%s]',donor_str);
if search_mode==1
    url = [service request probes donors];
elseif search_mode==2
    url = [service request donors structures probes];
end
if disp_url==1
disp(url);
end
%url='http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression%5Bdonors$eq12876%5D%5Bstructures$eq-----4006-----%5D%5Bprobes$eq1059351,1059352,1059353,1058915,1058914,1058916,1029156,1029155,1029150,1029149,1029146,1029142,1029138,1029136,1029135,1029134,1029133,1029126,1029125,1029124,1029123,1029122,1029121,1025133,1023492,1029129,1029128,1029127,1028516,1028515,1028514,1028513,1028512,1028511,1028510,1028509,1028508,1028507,1028506,1028505,1028504,1028503,1028502,1028484,1028483,1028478,1028477,1028474,1028473,1028472,1028471,1028470,1028469,1028466,1028465,1028464,1028462,1028461,1028458,1028457,1028456,1028455,1028454,1028453,1028452,1028451,1028450,1028448,1028447,1028446,1028445,1028444,1028442,1028441,1028440,1028439,1028438,1028437,1028434,1028433,1028432,1028431,1028430,1028429,1028428,1028427,1028426,1028425,1028424,1015329,1028436,1028435,1028501,1028500,1028499,1028498,1028497,1028494,1028493,1028492,1028491,1028490,1028489,1028488,1028487,1028486,1028485,1028548,1028547,1028542,1028536,1028535,1028534,1028527,1028526,1028518,1028517,1057976,1057975,1057974,1059659,1057967,1057966,1057965,1019360,1019310,1013314,1011485,1057962,1057961,1057960,1028346,1028345,1028344,1056548,1056547,1056546,1056545,1056544,1013593,1013489,1013359,1013692,1012504,1012098,1011320,1011137,1010770,1010725,1010443,1010639,1056550,1048801,1048800,1024804,1024803,1024812,1055392,1055381,1055378,1055374,1055373,1055372,1055361,1055357,1055356,1055355,1055354,1055347,1055346,1055345,1055342,1055377,1024617,1024606,1024594,1024593,1024592,1021988,1021983,1021982,1021981,1021979,1021975,1021974,1021973,1021972,1021970,1021969,1021968,1054369,1054368,1054364,1054362,1054360,1054359,1054358,1054357,1054353,1054352,1054350,1054349,1054348,1054347,1054345,1054344,1054343,1054342,1054341,1054340,1054339,1054338,1054337,1054336,1054335,1054334,1054380,1054379,1054375,1054372,1054370,1053289,1053288,1053287,1063261,1023147,1023146,1051465,1051464,1051461,1051460,1051459,1051458,1051457,1051489,1051483,1051472,1051466,1025754,1025749,1025748,1025747,1051004,1029261,1029258,1029257,1029256,1029255,1029254,1029253,1029252,1029251,1029250,1029249,1028377,1050663,1050662,1050661,1035832,1035831,1028612,1028611,1028610,1028609,1028608,1028607,1028606,1028605,1028604,1028603,1028602,1028601,1028634,1028633,1028631,1028627,1028626,1028625,1028620,1028619,1028618,1028617,1028615,1028614,1028613%5D'
status=0;
while status<1
    try
[str,status] = urlread(url);
    catch
        disp('lost connection. I start a new attempt!');
    end
end


json = parse_json(str);

% Pull out the samples, probes, and expression levels from the results.
samples = json.msg.samples;
nsamples = numel(samples);
probes = json.msg.probes;
explevels = zeros(nsamples,nprobes);
zscores = zeros(nsamples,nprobes);
for i=1:nprobes
    explevels(:,i) = cellfun(@str2num, probes{i}.expression_level);
    zscores(:,i) = cellfun(@str2num, probes{i}.z_score);
end        
