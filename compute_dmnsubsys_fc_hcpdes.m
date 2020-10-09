% This scripts computes average within- and between-connectivity for DMN
% susbsystem from DMN matrices derived from Schaefer 300 atlas timeseries
% of the HCP-DES dataset 


addpath('/share/leanew1/PANLab_Datasets/CONNECTOME/Leo/CIFTIMatlabReaderWriter');
addpath('/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/scripts');
addpath(genpath('/share/leanew1/PANLab_Datasets/CONNECTOME/Leo/gifti-1.6'))

setenv('LD_LIBRARY_PATH','/usr/lib')

clear

% Create DMN subsystem scores for HCP-DES
subs={'sub-CONN008', 'sub-CONN009', 'sub-CONN010', 'sub-CONN011', 'sub-CONN012', 'sub-CONN013', 'sub-CONN014', 'sub-CONN015', 'sub-CONN016', 'sub-CONN017', 'sub-CONN018', 'sub-CONN019', 'sub-CONN020', 'sub-CONN021', 'sub-CONN022', 'sub-CONN023', 'sub-CONN024', 'sub-CONN025', 'sub-CONN026', 'sub-CONN027', 'sub-CONN028', 'sub-CONN029', 'sub-CONN030', 'sub-CONN031', 'sub-CONN032', 'sub-CONN033', 'sub-CONN034', 'sub-CONN035', 'sub-CONN036', 'sub-CONN037', 'sub-CONN038', 'sub-CONN039', 'sub-CONN040', 'sub-CONN041', 'sub-CONN042', 'sub-CONN043', 'sub-CONN044', 'sub-CONN045', 'sub-CONN046', 'sub-CONN047', 'sub-CONN048', 'sub-CONN049', 'sub-CONN050', 'sub-CONN051', 'sub-CONN052', 'sub-CONN053', 'sub-CONN054', 'sub-CONN055', 'sub-CONN056', 'sub-CONN057', 'sub-CONN058', 'sub-CONN059', 'sub-CONN060', 'sub-CONN061', 'sub-CONN062', 'sub-CONN063', 'sub-CONN064', 'sub-CONN065', 'sub-CONN066', 'sub-CONN067', 'sub-CONN068', 'sub-CONN069', 'sub-CONN070', 'sub-CONN071', 'sub-CONN072', 'sub-CONN073', 'sub-CONN074', 'sub-CONN075', 'sub-CONN076', 'sub-CONN077', 'sub-CONN078', 'sub-CONN079', 'sub-CONN080', 'sub-CONN082', 'sub-CONN083', 'sub-CONN084', 'sub-CONN102', 'sub-CONN103', 'sub-CONN104', 'sub-CONN105', 'sub-CONN106', 'sub-CONN107', 'sub-CONN108', 'sub-CONN109', 'sub-CONN110', 'sub-CONN111', 'sub-CONN112', 'sub-CONN113', 'sub-CONN114', 'sub-CONN115', 'sub-CONN116', 'sub-CONN117', 'sub-CONN118', 'sub-CONN119', 'sub-CONN120', 'sub-CONN121', 'sub-CONN122', 'sub-CONN123', 'sub-CONN124', 'sub-CONN125', 'sub-CONN126', 'sub-CONN128', 'sub-CONN129', 'sub-CONN130', 'sub-CONN131', 'sub-CONN132', 'sub-CONN133', 'sub-CONN134', 'sub-CONN135', 'sub-CONN136', 'sub-CONN137', 'sub-CONN138', 'sub-CONN139', 'sub-CONN140', 'sub-CONN141', 'sub-CONN142', 'sub-CONN143', 'sub-CONN144', 'sub-CONN145', 'sub-CONN146', 'sub-CONN147', 'sub-CONN148', 'sub-CONN149', 'sub-CONN150', 'sub-CONN151', 'sub-CONN152', 'sub-CONN153', 'sub-CONN154', 'sub-CONN155', 'sub-CONN156', 'sub-CONN157', 'sub-CONN158', 'sub-CONN159', 'sub-CONN160', 'sub-CONN161', 'sub-CONN162', 'sub-CONN163', 'sub-CONN164', 'sub-CONN165', 'sub-CONN166', 'sub-CONN167', 'sub-CONN168', 'sub-CONN169', 'sub-CONN170', 'sub-CONN171', 'sub-CONN172', 'sub-CONN173', 'sub-CONN174', 'sub-CONN175', 'sub-CONN176', 'sub-CONN177', 'sub-CONN178', 'sub-CONN179', 'sub-CONN180', 'sub-CONN181', 'sub-CONN182', 'sub-CONN183', 'sub-CONN184', 'sub-CONN185', 'sub-CONN186', 'sub-CONN187', 'sub-CONN188', 'sub-CONN189', 'sub-CONN190', 'sub-CONN191', 'sub-CONN192', 'sub-CONN193', 'sub-CONN194', 'sub-CONN195', 'sub-CONN196', 'sub-CONN197', 'sub-CONN198', 'sub-CONN199', 'sub-CONN200', 'sub-CONN201', 'sub-CONN202', 'sub-CONN203', 'sub-CONN204', 'sub-CONN205', 'sub-CONN206', 'sub-CONN207', 'sub-CONN208', 'sub-CONN209', 'sub-CONN210', 'sub-CONN211', 'sub-CONN212', 'sub-CONN213', 'sub-CONN214', 'sub-CONN215', 'sub-CONN216', 'sub-CONN217', 'sub-CONN218', 'sub-CONN219', 'sub-CONN220', 'sub-CONN221', 'sub-CONN222', 'sub-CONN223', 'sub-CONN224', 'sub-CONN225', 'sub-CONN226', 'sub-CONN227', 'sub-CONN228', 'sub-CONN230', 'sub-CONN231', 'sub-CONN232', 'sub-CONN234', 'sub-CONN235', 'sub-CONN236', 'sub-CONN237'};

% Directories
ciftidir='/oak/stanford/groups/leanew1/ltozzi/hcpdes/aim2_project/CIFTI_data/';
dmnsubsys_outputdir='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/dmn_subsystems/';
matrices_outputdir='/oak/stanford/groups/leanew1/ltozzi/hcpdes/dmn_subsystems_project/matrices/';

% Indeces for DMN subsystems in the schaefer 300 atlas
coreidx=[111 112 113 114 115 116 117 118 119 120 121 271 272 273 274 275 276 277 278 279 280 281 282];
dmpfcidx=[122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 283 284 285 286 287 288 289];
mtlidx=[140 141 142 143 144 145 290 291 292 293];

dmnidx=[coreidx dmpfcidx mtlidx];

incomplete=zeros(length(subs), 1);
fd_vols_all=nan(length(subs), 1);
allmats=struct();

fdthresh=0.25;

for sub=1:length(subs)
    
    disp(subs{sub})
    
    % check if data exists
    run1=strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-00_task-rest_acq-mb_dir-pe0_space-fsLR_den-32k_bold_AROMAclean_sch300N17.ptseries.nii');
    run2=strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-00_task-rest_acq-mb_dir-pe1_space-fsLR_den-32k_bold_AROMAclean_sch300N17.ptseries.nii');
    run3=strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-01_task-rest_acq-mb_dir-pe0_space-fsLR_den-32k_bold_AROMAclean_sch300N17.ptseries.nii');
    run4=strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-01_task-rest_acq-mb_dir-pe1_space-fsLR_den-32k_bold_AROMAclean_sch300N17.ptseries.nii');
    
    if ~(exist(run1) && exist(run2) && exist(run3) && exist(run4))
        
        disp('Missing data')
        incomplete(sub)=1;
        
    else
        % load each run, filter and concatenate
        cif=ciftiopen(run1, 'wb_command');
        run1_data=highpass(cif.cdata', 0.01, 1/0.710);
        cif=ciftiopen(run2, 'wb_command');
        run2_data=highpass(cif.cdata', 0.01, 1/0.710);
        cif=ciftiopen(run3, 'wb_command');
        run3_data=highpass(cif.cdata', 0.01, 1/0.710);
        cif=ciftiopen(run4, 'wb_command');
        run4_data=highpass(cif.cdata', 0.01, 1/0.710);
        
        rest_data=[run1_data;run2_data;run3_data;run4_data];
        
        % load fd for each run
        reg=readtable(strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-00_task-rest_acq-mb_dir-pe0_desc-confounds_regressors.tsv'), 'FileType', 'text');
        fd1=str2double(reg.framewise_displacement);
        fd1(1)=0;
        reg=readtable(strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-00_task-rest_acq-mb_dir-pe1_desc-confounds_regressors.tsv'), 'FileType', 'text');
        fd2=str2double(reg.framewise_displacement);
        fd2(1)=0;
        reg=readtable(strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-01_task-rest_acq-mb_dir-pe0_desc-confounds_regressors.tsv'), 'FileType', 'text');
        fd3=str2double(reg.framewise_displacement);
        fd3(1)=0;
        reg=readtable(strcat(ciftidir, subs{sub}, '/func/', subs{sub}, '_ses-01_task-rest_acq-mb_dir-pe1_desc-confounds_regressors.tsv'), 'FileType', 'text');
        fd4=str2double(reg.framewise_displacement);
        fd4(1)=0;
        
        rest_fd=[fd1;fd2;fd3;fd4];
        
        % censor volumes
        if sum(rest_fd>fdthresh)/size(rest_fd, 1)>0.25
            disp(strcat(subs{sub}, {' '}, 'too much motion'))
            incomplete(sub)=1;
        else
            
            % scrub if FD>fdthresh
            rest_data(rest_fd>fdthresh,:)=[];
            fd_vols_all(sub)=sum(rest_fd>fdthresh);
            
            % save connectivity matrix
            allmats(sub).id=subs{sub};
            allmats(sub).dmn=atanh(corr(rest_data(:, dmnidx)));
            allmats(sub).nfd=sum(rest_fd>fdthresh);
            
            % calculate connectivity within each subsystem
            mat_core=atanh(corr(rest_data(:, coreidx)));
            corevec(sub)=mean(unpackconnmat(mat_core));
            mat_dmpfc=atanh(corr(rest_data(:, dmpfcidx)));
            dmpfcvec(sub)=mean(unpackconnmat(mat_dmpfc));
            mat_mtl=atanh(corr(rest_data(:, mtlidx)));
            mtlvec(sub)=mean(unpackconnmat(mat_mtl));
            
            % calculate connectivity between each subsystem
            core_dmpfcvec(sub)=mean(mean(atanh(corr(rest_data(:, coreidx), rest_data(:, dmpfcidx)))));
            core_mtlvec(sub)=mean(mean(atanh(corr(rest_data(:, coreidx), rest_data(:, mtlidx)))));
            dmpfc_mtlvec(sub)=mean(mean(atanh(corr(rest_data(:, dmpfcidx), rest_data(:, mtlidx)))));
        end
    end
end


% create table

subs(logical(incomplete))=[];
fd_vols_all(logical(incomplete))=[];
corevec(logical(incomplete))=[];
dmpfcvec(logical(incomplete))=[];
mtlvec(logical(incomplete))=[];
core_dmpfcvec(logical(incomplete))=[];
core_mtlvec(logical(incomplete))=[];
dmpfc_mtlvec(logical(incomplete))=[];

allmats(logical(incomplete))=[];

nwstable_hcpdes=table();
nwstable_hcpdes.ID = subs';
nwstable_hcpdes.fdvols = fd_vols_all;
nwstable_hcpdes = addvars(nwstable_hcpdes, corevec','NewVariableNames', 'Core');
nwstable_hcpdes = addvars(nwstable_hcpdes, dmpfcvec','NewVariableNames', 'DMPFC');
nwstable_hcpdes = addvars(nwstable_hcpdes, mtlvec','NewVariableNames', 'MTL');
nwstable_hcpdes = addvars(nwstable_hcpdes, core_dmpfcvec','NewVariableNames', 'Core_DMPFC');
nwstable_hcpdes = addvars(nwstable_hcpdes, core_mtlvec','NewVariableNames', 'Core_MTL');
nwstable_hcpdes = addvars(nwstable_hcpdes, dmpfc_mtlvec','NewVariableNames', 'DMPFC_MTL');

% save table
writetable(nwstable_hcpdes, strcat(dmnsubsys_outputdir, 'DMNsubsys_fc_hcpdes.csv'));

% save structure
save(strcat(matrices_outputdir, 'connmatstack_dmn_hcpdes.mat'), 'allmats');
