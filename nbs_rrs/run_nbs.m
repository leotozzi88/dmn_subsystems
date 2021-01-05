% Prep NBS

clear

% Load the clinical info 
hcpdestab=readtable('../../clinical/hcpdes.csv');

% Load the matrix stack
load('../../matrices/connmatstack_dmn_hcpdes.mat');

% Create NBS matrix stack for subs who have a matrix and rrs
i=1;
for sub=1:size(allmats, 2)
    subid = getfield(allmats, {sub},'id');
    if ismember(subid, hcpdestab.ID) && ~isempty(hcpdestab(strcmp(hcpdestab.ID, subid), 'rrs_total'))
        mat=getfield(allmats, {sub},'dmn');
        if ~isnan(mean(unpackconnmat(mat))) % double check for nans
            finalsubslist{i}=subid;
            matstack(:,:, i) = getfield(allmats, {sub},'dmn');
            i=i+1;
        end
    end
end

save('dmn_connmats_hcpdes_rrs', 'matstack')

% Create design matrix 
hcpdestab(~ismember(hcpdestab.ID, finalsubslist), :)=[];

desmat=hcpdestab(:, {'bio_sex','demo_age', 'rrs_total'});
writetable(desmat, 'hcpdes_nbs_rrs_total_desmat.txt', 'Delimiter', '\t', 'WriteVariableNames', 0);

desmat=hcpdestab(:, {'bio_sex','demo_age', 'reflection_total'});
writetable(desmat, 'hcpdes_nbs_rrs_reflection_desmat.txt', 'Delimiter', '\t', 'WriteVariableNames', 0);

desmat=hcpdestab(:, {'bio_sex','demo_age', 'brooding_total'});
writetable(desmat, 'hcpdes_nbs_rrs_brooding_desmat.txt', 'Delimiter', '\t', 'WriteVariableNames', 0);

desmat=hcpdestab(:, {'bio_sex','demo_age', 'deprelated_total'});
writetable(desmat, 'hcpdes_nbs_rrs_deprelated_desmat.txt', 'Delimiter', '\t', 'WriteVariableNames', 0);


% Run NBS

clear nbs_results

tidx=1;
for thresh=1:0.2:5
    
    UI=struct();
    UI.method.ui='Run NBS';
    UI.test.ui='t-test';
    UI.size.ui='Intensity';
    UI.thresh.ui=num2str(thresh);
    UI.perms.ui='10000';
    UI.alpha.ui='0.05';
    UI.contrast.ui='[0 0 -1]';
    UI.design.ui='hcpdes_nbs_rrs_reflection_desmat.txt';
    UI.exchange.ui='';
    UI.matrices.ui='dmn_connmats_hcpdes_rrs_cn.mat';
    UI.node_coor.ui='';
    UI.node_label.ui='';
    
    NBSrun(UI, '');
    global nbs
    
    % store the nbs result
    nbs_results{tidx}=nbs;
    tidx=tidx+1;
    
end

% Plot results
clear nnws
for n=1:length(nbs_results)

    nnws(n)=nbs_results{1, n}.NBS.n;
    
end

plot(nnws)
ylim([-1 2])
