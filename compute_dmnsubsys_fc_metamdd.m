% This scripts computes average within- and between-connectivity for DMN
% susbsystem from DMN matrices derived from Power atlas timeseries of the Meta-MDD
% dataset 

%%% CALCULATE AVERAGE CONNECTIVITY PER SUBSYSTEM

clear all

% Directories
outputdir='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/dmn_subsystems/';
inputdir='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/matrices/';

load(strcat(inputdir, 'connmatstack_dmn_metamdd.mat'))
load(strcat(inputdir, 'group_dmn_metamdd.mat'))
load(strcat(inputdir, 'subslist_dmn_metamdd.mat'))
load(strcat(inputdir, 'fd_volumes_dmn_metamdd.mat'))

% create table
Gnew=G;
Gnew(Gnew==2)=-1;
nwstable=table();
nwstable.ID = subs';
nwstable.Group=Gnew;
nwstable.fdvols = fd_vols_all;

% Subnetwork assignment for sub-systems within the DMN
P_dmnsubsys=[3 0 1 0 0 2 3 2 2 2 1 1 1 1 3 1 1 0 0 3 1 1 2 2 2 1 2 2 1 1 2 1 1 1 1 1 0 1 1 2 1 2 2 0 2 2 1 2 3 3 3 0 2 2 0 2 2 2];
dmnsubsys_names={'Core', 'DMPFC', 'MTL'};

% Add within sub-systems connectivity to table

for nw=1:max(P_dmnsubsys)
    
    nw_vec=nan(length(subs), 1);
    
    for sub=1:size(stack, 3)
        
        submat=stack(:, :, sub);
        nwmat=submat(P_dmnsubsys==nw, P_dmnsubsys==nw);
        mask=logical(triu(ones(size(nwmat)), 1));
        nw_vec(sub)=mean(nwmat(mask));
        
    end
    
    % update table
    nwstable = addvars(nwstable, nw_vec,'NewVariableNames', dmnsubsys_names{nw});
    
end

% save table
writetable(nwstable, strcat(outputdir, 'within_DMNsubsys_fc_metamdd.csv'));


%%% CALCULATE AVERAGE CONNECTIVITY BETWEEN NETWORK


% create table
btwnwstable=table();
btwnwstable.ID = subs';
btwnwstable.Group=Gnew;
btwnwstable.fdvols = fd_vols_all;

for nw1=1:max(P_dmnsubsys)
    
    for nw2=1:max(P_dmnsubsys)
        
        nw_vec=nan(length(subs), 1);
        
        if nw1~=nw2
            
            for sub=1:size(stack, 3)
                
                submat=stack(:,:,sub);
                nwmat=submat(P_dmnsubsys==nw1, P_dmnsubsys==nw2);
                nw_vec(sub)=mean(nwmat(:));
                
            end
            
            % check if a duplicate column exists
            if ~ismember(strcat(dmnsubsys_names{nw2}, '-', dmnsubsys_names{nw1}), btwnwstable.Properties.VariableNames)
            
            % update table
            btwnwstable = addvars(btwnwstable, nw_vec,'NewVariableNames', strcat(dmnsubsys_names{nw1}, '_', dmnsubsys_names{nw2}));
            end
        end
    end
end

% save table
writetable(btwnwstable, strcat(outputdir, 'between_DMNsubsys_fc_metamdd.csv'));