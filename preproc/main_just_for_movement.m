%% This script will be used for the first time screening of the data (only give you movement info)
%This script was created by Professor Baxter Rogeres (VUIIS), but is heavily modified based on our lab pipeline by Jin Wang 3/7/2019
%(1) realignment to mean, reslice the mean.
%(2) segment anatomical image to TPM template. We get a deformation file "y_filename" and this is used in normalisation step to normalize all the
%    functional data and the mean functional data.
%(3) Then we make a skull-striped anatomical T1 (based on segmentation) and coregister mean functional data (and all other functional data) to the anatomical T1.
%(4) Smoothing.
%(5) Art_global. It calls the realignmentfile (the rp_*.txt) to do the interpolation. This step identifies the bad volumes(by setting scan-to-scan movement
%    mv_thresh =1.5mm and global signal intensity deviation Percent_thresh= 4 percent, any volumes movement to reference volume, which is the mean, >5mm) and repair
%    them with interpolation. This step uses art-repair art_global.m function (the subfunctions within it are art_repairvol, which does repairment, and art_climvmnt, which identifies volumes movment to reference.
%(6) We use check_reg.m to see how well the meanfunctional data was normalized to template by visual check.
%(7) collapse all the 3d files into 4d files for saving space. You can decide whether you want to delete the product of not later.
%Before you run this script, you should make sure that your data structure is as expected.



global CCN;
addpath(genpath('/dors/booth/JBooth-Lab/BDL/jinwang/SemPhon_7_8/scripts')); %This is the code path
spm_path='/dors/booth/JBooth-Lab/BDL/LabCode/typical_data_analysis/spm12_elp'; %This is your spm path
%tpm='/dors/gpc/JamesBooth/JBooth-Lab/BDL/LabCode/typical_data_analysis/templates_cerebroMatic/mw_com_prior_Age_0081.nii'; %This is your template path
addpath(genpath(spm_path));
root='/dors/booth/JBooth-Lab/BDL/jinwang/SemPhon_7_8'; %This is your project folder
test_subjects={'sub-5003'};
run_script=1; % 1 is to run test_subjects, 2 is to run all the rest of the subjects in the preprocessed folder that's not specified in test_subjects. 

%if user did not fully specify the subjects list, then it will read in all the
%data except the ones specified in the subjects list
listing=dir([root '/' CCN.preprocessed]);
all_list=extractfield(listing,'name');
index=strfind(all_list,'sub');
idx=find(not(cellfun('isempty',index)));
subjects_all=all_list(idx);
if run_script ==1
    subjects=test_subjects;
else
    subjects=subjects_all(~ismember(subjects_all,test_subjects));
end

CCN.preprocessed_folder='preprocessed'; %This is your data folder needs to be preprocessed
CCN.session='ses-7'; %This is the time point you want to analyze. If you have two time points, do the preprocessing one time point at a time by changing ses-T1 to ses-T2.
CCN.func_folder='sub*'; % This is your functional folder name
CCN.func_pattern='sub*.nii'; %This is your functional data name
%CCN.anat_pattern='ses-T1_T1w*.nii'; %This is your anat data name


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%shouldn't be modified below%%%%%%%%%%%%%%%
% Initialize
%addpath(spm_path);
spm('defaults','fmri');
spm_jobman('initcfg');
spm_figure('Create','Graphics','Graphics');

% Dependency and sanity checks
if verLessThan('matlab','R2013a')
    error('Matlab version is %s but R2013a or higher is required',version)
end

req_spm_ver = 'SPM12 (6225)';
spm_ver = spm('version');
if ~strcmp( spm_ver,req_spm_ver )
    error('SPM version is %s but %s is required',spm_ver,req_spm_ver)
end
try
    %Start to preprocess data from here
    for i=1:length(subjects)
        fprintf('work on subject %s\n', subjects{i});
        CCN.subj_folder=[root '/' CCN.preprocessed_folder '/' subjects{i}];
        out_path=[CCN.subj_folder '/' CCN.session];
        CCN.func_f='[subj_folder]/[session]/func/[func_folder]/';
        func_f=expand_path(CCN.func_f);
        func_file=[];
        for m=1:length(func_f)
            func_file{m}=expand_path([func_f{m} '[func_pattern]']);
        end
        
        
        Percent_thresh= 4; %global signal intensity change
        mv_thresh =1.5; % scan-to-scan movement
        MVMTTHRESHOLD=5; % movement to reference,see in art_clipmvmt
        
        %expand 4d functional data to 3d data
        func_vols=cell(length(func_file),1);
        for x=1:length(func_file)
            hdr = load_nii_hdr( char(func_file{x}) );
            nscan = hdr.dime.dim(5);
            expand_nii_scan(char(func_file{x}),1:nscan)
            [func_p,func_n] = fileparts(char(func_file{x}));
            func_vols{x}=cellstr(spm_select('ExtFPList',func_p, [func_n '_0.*\.nii$']));
            [rfunc_file, rp_file] = realignment_byrun(func_vols{x}, out_path);%this will run the realignment by run by run
            for jj=1:length(rfunc_file)
                swfunc_list{jj}=erase(rfunc_file{jj},',1');
            end
            art_global_jin(char(swfunc_list),rp_file,4,1,Percent_thresh,mv_thresh,MVMTTHRESHOLD);
            for ii=1:length(swfunc_list)
                delete(swfunc_list{ii});
                [p,n,ext]=fileparts(swfunc_list{ii});
                delete([p '/v' n ext]);
            end
            [p,n,ext]=fileparts(swfunc_list{1});
            delete([p '/mean' n ext]);
            delete([p '/rp_' n '.txt']);
        end
    end
catch e
    rethrow(e)
    %display the errors
end
