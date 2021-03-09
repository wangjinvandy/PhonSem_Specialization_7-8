%% First level analysis, written by Jin Wang 3/15/2019
% You should define your conditions, onsets, duration, TR.
% The repaired images will be deweighted from 1 to 0.01 in the first level
% estimation (we changed the art_redo.m, which uses art_deweight.txt as default to deweight the scans to art_redo_jin.txt, which we uses the art_repaired.txt to deweight scans).
% The difference between art_deiweghted.txt and art_repaired.txt is to the
% former one is more wide spread. It not only mark the scans which was repaired to be deweighted, but also scans around it to be deweighted.
% The 6 movement parameters we got from realignment is added into the model regressors to remove the small motion effects on data.

% Make sure you run clear all before running this code. This is to clear
% all existing data structure which might be left by previous analysis in
% the work space.

%This code now is only for bids project that has unique runs (not like ELP project that we sometimes have multiple runs for one task run).

addpath(genpath('/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/SemPhon_7_8/scripts/3firstlevel')); % the path of your scripts
spm_path='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/LabCode/typical_data_analysis/spm12_elp'; %the path of spm
addpath(genpath(spm_path));

%define your data path
data=struct();
root='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/SemPhon_7_8';  %your project path
subjects={};
data_info='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/jinwang/SemPhon_7_8/SemPhon_7_8.xlsx'; %final_sample sheet, the subject column should have sub plus subject number.(e.g. sub-5003)
if isempty(subjects)
    M=readtable(data_info);
    subjects=M.subjects;
end

analysis_folder='analysis'; % the name of your first level modeling folder
model_deweight='deweight'; % the deweigthed modeling folder, it will be inside of your analysis folder

global CCN
CCN.preprocessed='preprocessed'; % your data folder
CCN.session='ses*'; % the time points you want to analyze
CCN.func_pattern='sub*'; % the name of your functional folders
CCN.file='vs6_wsub*bold.nii'; % the name of your preprocessed data (4d)
CCN.files='vs6_wsub*bold_0*'; % the name of your preprocessed data (3d) expanded one.
CCN.rpfile='rp_*.txt'; %the movement files
events_file_exist=0; % 1 means you did not clean your events.tsv, 0 means you cleaned the events.tsv in your preprocessed folder
bids_folder='/gpfs51/dors2/gpc/JamesBooth/JBooth-Lab/BDL/ELP/bids'; % if you assign 0 to events_file_exist, then you mask fill in this path, so it can read events.tsv file for individual onsets from bids folder

%%define your task conditions, be sure it follows the sequence of the
%%output from count_repaired_acc_rt.m.
conditions=[];
conditions{1}={'P_O' 'P_R' 'P_SC' 'P_U'};
conditions{2}={'P_O' 'P_R' 'P_SC' 'P_U'};
conditions{3}={'S_H' 'S_L' 'S_SC' 'S_U'};
conditions{4}={'S_H' 'S_L' 'S_SC' 'S_U'};

%duration
dur=0; %I think all projects in BDL are event-related, so I hard coded the duration as 0.

%TR
TR=1.25; %ELP project

%define your contrasts, make sure your contrasts and your weights should be
%matched.
contrasts={'onsetrhyme_vs_perceptual' ...
    'weakstrong_vs_perceptual'...
    'onsetrhyme_vs_perc_VS_weakstrong_vs_perc' ...
    'weakstrong_vs_perc_VS_onsetrhyme_vs_perc' ...
    'onset_vs_rhyme' ...
    'weak_vs_strong'};
onsetrhyme_vs_perceptual=[1 1 -2 0];
weakstrong_vs_perceptual=[1 1 -2 0];
onset_vs_rhyme=[1 -1 0 0];
weak_vs_strong=[-1 1 0 0];

%adjust the contrast by adding six 0s into the end of each session
rp_w=zeros(1,6);
empty=zeros(1,10);
weights={[onsetrhyme_vs_perceptual rp_w onsetrhyme_vs_perceptual rp_w empty empty]...
    [empty empty weakstrong_vs_perceptual rp_w weakstrong_vs_perceptual rp_w] ...
    [onsetrhyme_vs_perceptual rp_w onsetrhyme_vs_perceptual rp_w -1*(weakstrong_vs_perceptual) rp_w -1*(weakstrong_vs_perceptual) rp_w] ...
    [-1*(onsetrhyme_vs_perceptual) rp_w -1*(onsetrhyme_vs_perceptual) rp_w weakstrong_vs_perceptual rp_w weakstrong_vs_perceptual rp_w] ...
    [onset_vs_rhyme rp_w onset_vs_rhyme rp_w empty empty] ...
    [empty empty weak_vs_strong rp_w weak_vs_strong rp_w]};

%%%%%%%%%%%%%%%%%%%%%%%%Do not edit below here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check if you define your contrasts in a correct way
if length(weights)~=length(contrasts)
    error('the contrasts and the weights are not matched');
end

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

%Start to analyze the data from here
try
    for i=1:length(subjects)
        fprintf('work on subject %s', subjects{i});
        CCN.subject=[root '/' CCN.preprocessed '/' subjects{i}];
        %specify the outpath,create one if it does not exist
        out_path=[CCN.subject '/' analysis_folder];
        if ~exist(out_path)
            mkdir(out_path)
        end
         
        %specify the deweighting spm folder, create one if it does not exist
        model_deweight_path=[out_path '/' model_deweight];
        if exist(model_deweight_path,'dir')~=7
            mkdir(model_deweight_path)
        end
        
        %find folders in func
        CCN.functional_dirs='[subject]/[session]/func/[func_pattern]/';
        functional_dirs=expand_path(CCN.functional_dirs);
          
        %load 6 movement parameters
        mv=[];
        rp_file=expand_path([CCN.functional_dirs '[rpfile]']);
        for i=1:length(rp_file)
            rp=load(rp_file{i});
            mv{i}=rp;
        end
        data.mv=mv;
        
        %load the functional data and event onsets
        swfunc=[];
        P=[];
        onsets=[];
        for j=1:length(functional_dirs)
            %load the functional data
            P{j}=char(expand_path([functional_dirs{j} '[files]']));
            if isempty(P{j}) % if there is no expanded data, then expand the preprocessed data
                file=expand_path([functional_dirs{j} '[file]']);
                hdr = load_nii_hdr(char(file));
                nscan = hdr.dime.dim(5);
                expand_nii_scan(char(file),1:nscan);
                P{j}=char(expand_path([functional_dirs{j} '[files]']));
            end
            [m,n]=size(P{j});
            for ii=1:m
                swfunc{j}(ii,:)=[P{j}(ii,:) ',1'];
            end
            %load the event onsets
            if events_file_exist==1
                [p,event_n]=fileparts(functional_dirs{j}(1:end-1));
                event_file=[functional_dirs{j} event_n(1:end-4) 'events.tsv'];
            elseif events_file_exist==0
                [p,event_n]=fileparts(functional_dirs{j}(1:end-1));
                [q,session]=fileparts(fileparts(p));
                [~,this_subject]=fileparts(q);
                event_file=[bids_folder '/' this_subject '/' session '/func/' event_n(1:end-4) 'events.tsv'];
            end
            event_data=tdfread(event_file);
            cond=unique(event_data.trial_type,'row');
            for k=1:size(cond,1)
            onsets{j}{k}=event_data.onset(sum((event_data.trial_type==cond(k,:))')'==4);
            end
            
        end
        data.swfunc=swfunc;
        
        %pass the experimental design information to data
        data.conditions=conditions;
        data.onsets=onsets;
        data.dur=dur;
        
        %run the firstlevel modeling and estimation (with deweighting)
        mat=firstlevel(data, out_path, TR, model_deweight_path);
        origmat=[out_path '/SPM.mat'];
        %run the contrasts
        contrast(origmat,contrasts,weights);
        contrast(mat,contrasts,weights);
        
%         %delete the expanded vs6_wtask*.nii data in your preprocessing folder
%         for de=1:length(swfunc)
%                     [q,p]=size(swfunc{de});
%             for del=1:q
%                 file_to_delete=swfunc{de}(del, 1:end-2);
%                 delete(file_to_delete);
%             end
%         end
        
    end
    
catch e
    rethrow(e)
    %display the errors
end