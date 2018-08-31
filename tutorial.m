%{

%% Required software

in this tarball you will find the 2 sessions from the same subject (i.e. 101295 is session 1 and
101363 is session 2).

The two sessions are aligned, and the structural of the 1st session is also
aligned within a subject.

The 18 first images (after applying reorderBrainBetas are
the subject's own images from his photo-album, the following 18 images are
the other paired subject's own images (we will not analyse that subject's fMRI data today
and then the last 36 images are the general images that every subject from every
pair saw.

the similarity judgments RDM are contained in the
pair1_subj1_extra.mat in the structure behav_and_icons.RDM

Here are some more notes on the stimulus indices:
 
OwnBodyParts = 1:3; OwnFaces = 4:8; OwnPet = 9; OwnPlaces = 10:15; OwnObjects=16:18;
OtherBodyParts = 19:21; OtherFaces = 22:26; OtherPet = 27; OtherPlaces = 28:33; OtherObjects=34:36;
GeneralBodyParts = 37:44; GeneralFaces = 45:52; GeneralPets = [53 54]; GeneralPlaces=55:66; GeneralObjects=67:72;

note that own refer to the subjects' own images, and other to the paired
subject's own images, and general to objects that were seen by all subjects.

bodies = [OwnBodyParts OtherBodyParts GeneralBodyParts];
nbodies = length(bodies);
faces = [OwnFaces OwnPet OtherFaces OtherPet GeneralFaces GeneralPets];
nfaces = length(faces);
places = [OwnPlaces OtherPlaces GeneralPlaces];
nplaces = length(places);
objects = [OwnObjects OtherObjects GeneralObjects];
nobjects = length(objects);

animates = [faces bodies];
inanimates = [places objects];


% ALREADY DONE
% 1 - preprocess and split the data
% 2 - estimate single-subject activity patterns

% LINEAR SVM
% 2 - review activity-pattern estimation
        inspect the design matrices for the runs in both sessions
        load the single-subject activity patterns (t-patterns)
% 3 - select voxels
        constrain the t-patterns to the hIT masks (actually called  VOT)
% 4 - train the classifier
%       train a linear SVM classifier to distinguish animate from inanimate
%       objects
% 5 - test the classifier
% 6 - statistical inference
%       use a condition-label permutation test to determine whether the
%       classifier performs above chance in this particular subject

% RSA
% 7 - generate an RDM for each session
% 8 - compare RDMs across sessions
% 9 - relate the session-averaged hIT RDM to the similarity judgments RDM

%}

clear all
clc
[workingdir,file,ext] = fileparts(which('rsa_tutorial'));

analyse.lSVM=0;
analyse.RSA=1;


%% required software --> SPM, GLMdenoise, libSVM, RSA toolbox
addpath(genpath(fullfile(workingdir,'software','GLMdenoise-1.4')));
addpath(genpath(fullfile(workingdir,'software','rsatoolbox')));
addpath(genpath(fullfile(workingdir,'software','spm8')));
addpath(fullfile(workingdir,'software','libsvm-mat-2.87-1'));


%% control parameters

nRuns        = 9;   % there was a total of 9 runs per subject and session
nConditions  = 72;  % there were 72 different object conditions in the experiment
nVolumes     = 216; % there was 216 EPI volumes acquired per run
tr           = 2;   % the EPI volumes were acquired every 2 seconds
stimdur      = 1;   % the stimuli were presented on screen for 1 second

% untar the pair 1 subject 1 tarball
subjecttar   = 'pair1_subj1.tar.gz';

if ~(exist('CBU101295','dir'))
    system(sprintf('tar -xvzf %s',subjecttar));
end

% the untarring will extract the two following folders
sessions     = {'CBU101295';'CBU101363'};
nSessions    = numel(sessions);

% load the subject's image ordering (the same ordering is valid for both sessions)
load('pair1_subj1_extra.mat','reorderBrainBetas')

% pre allocate cell arrays
tempOnsets   = cell(1,nRuns);
design       = cell(1,nRuns);
tpatterns        = cell(1,nSessions);

% the idendity of the subject will be session 1 (useful when loading mask)
thisSubject = sessions{1};

% load hIT mask
mask = logical(spm_read_vols(spm_vol(fullfile(workingdir,thisSubject,'masks','VOT','bilateral.VOT.nii'))));

% vectorize the mask
vmask = squish(mask,3);


% loop over the 2 sessions
for sessionI = 1:nSessions
    
    thisSess = sessions{sessionI};
    
    % verify that the nVolumes match
    % V=spm_vol(sprintf('%s/run1/4d.nii',thisSess));
    % assert(size(V,1)==nVolumes)
    
    logfile_dir = sprintf('%s/logfiles',thisSess);
    
    % prepare design matrices
    tempOnsets = cell(1,nRuns);
    tempDesign = cell(1,nRuns);
    
    for runI=1:nRuns
        
        fileName=fullfile(logfile_dir,[thisSess '.iRSAexperiment.' num2str(runI) '.txt']);
        %read the scanning log file
        [trialI,responseType_CrHMFa,rt,starttime,stimType,stimNumber,stimName] = textread(fileName, '%d %d %f %f %d %d %s','headerlines',7);
        
        starttime = (starttime/1000);
        % the first session of the 1st subject had a bug with the onsets,
        % which you can correct like follows:
        if strcmp(thisSess,'CBU101295')
            starttime = 32:4:412;
        else
            starttime = round(starttime);
        end
        
        for condI=1:nConditions
            tempOnsets{condI}=starttime(stimNumber==condI & stimType==0);
            % tempOnsets{condI}=starttime(stimNumber==condI & (stimType==0 | stimType ==1));
        end
        % anomalies
        tempOnsets{nConditions+1}=starttime(stimType==1);
        
        % reorder the onsets according to conceptual slots
        for condI=1:nConditions
            tempDesign{runI}{condI}=tempOnsets{reorderBrainBetas(condI)};
        end
        tempDesign{runI}{nConditions+1}=tempOnsets{nConditions+1}; % deal with the anomalous trials
    end
    design = cell(1,nRuns);
    for runI=1:nRuns
        design{runI} = zeros(nVolumes,nConditions+1);
        for condI=1:nConditions+1
            design{runI}(tempDesign{runI}{condI}/2+1,condI) = 1;
        end
    end
    
    % plot one run's design matrix
    figure(1);
    set(gcf,'Position',[100 100 400 800],'Color','w')
    imagesc(design{1});colormap(gray)
    ylabel('\bf{number of volumes}')
    xlabel('\bf{number of conditions}')
    title(sprintf('design matrix for run 1 session %d',sessionI),'FontSize',14,'Fontweight','bold')
    
    if ~(exist(fullfile(workingdir,'tpatterns','tpatterns.mat'),'file'))
        data = cell(1,nRuns);
        % load the nifti data;
        for runI=1:nRuns
            fprintf('***\t importing EPI time-series for run %d \t***\n',runI)
            thisRun = fullfile(workingdir,thisSess,sprintf('run%d',runI));
            data{runI} = single(spm_read_vols(spm_vol(fullfile(thisRun,'4d.nii'))));
        end
        
        % run GLMdenoise
        results = GLMdenoisedata(design,data,stimdur,tr, ...
            'optimize',[],struct('numboots',100,'numpcstotry',20,'wantparametric',1), ...
            []);
        
        % limit the betas to the valid conditions
        modelmd = results.modelmd{2}(:,:,:,1:nConditions);
        % limit the standard errors to the valid conditions
        modelse = results.modelse{2}(:,:,:,1:nConditions);
        % get the pooled error
        poolse  = sqrt(mean(modelse.^2,4));
        % normalise the betas by the pooled error to get t-patterns
        modelmd = bsxfun(@rdivide,modelmd,poolse);
        
        % show the unmasked patterns
        figure(1);
        set(gcf,'Position',[100 100 600 800],'Color','w')
        subplot(3,1,1)
        imagesc(makeimagestack(mean(modelmd,4),-1))
        title('\bf{unmasked mean t-pattern}')
        subplot(3,1,2)
        imagesc(makeimagestack(mask,-1))
        title('\bf{hIT mask}')
        subplot(3,1,3)
        imagesc(makeimagestack(mean(modelmd,4).*mask,-1))
        title('\bf{mean t-pattern masked to hIT}')
        
        vmodelmd     = squish(modelmd,3);
        
        % masked t-pattenrs will be nvoxels x conditions
        tpatterns{sessionI}    = vmodelmd(vmask,:);
    end
    
end

% save the data if not existing already otherwise load it
if ~(exist(fullfile(workingdir,'tpatterns','tpatterns.mat'),'file'))
    % save the t-patterns
    save(fullfile(workingdir,'tpatterns','tpatterns.mat'),'tpatterns')
else
    % if the t-patterns already exist --> load it
    load(fullfile(workingdir,'tpatterns','tpatterns.mat'))
end


%% animate vs. inanimate classification using linear SVM

if analyse.lSVM
    % define animate and inanimate index vectors
    OwnBodyParts = 1:3; OwnFaces = 4:8; OwnPet = 9; OwnPlaces = 10:15; OwnObjects=16:18;
    OtherBodyParts = 19:21; OtherFaces = 22:26; OtherPet = 27; OtherPlaces = 28:33; OtherObjects=34:36;
    GeneralBodyParts = 37:44; GeneralFaces = 45:52; GeneralPets = [53 54]; GeneralPlaces=55:66; GeneralObjects=67:72;
    
    animates=[OwnBodyParts OwnFaces OwnPet OtherBodyParts OtherFaces OtherPet GeneralBodyParts GeneralFaces GeneralPets];
    inanimates=[OwnPlaces OwnObjects OtherPlaces OtherObjects GeneralPlaces GeneralObjects];
    
    % control variables
    libSVMsettings='-s 1 -t 0'; % nu-SVM, linear
    nRandomisations=1000;
    %rmpath('/hpc-software/matlab/r2009a/toolbox/bioinfo/biolearning/'); % to make sure libSVM code is used (not strictly necessary: matlab svmtrain yields exactly same model)
    
    % linear SVM
    cvFolds=[1 2; 2 1]; % columns = folds, row 1 = session used for training, row 2 = session used for testing
    for foldI=1:size(cvFolds,2)
        % define training and test data sets
        tpatternsTrain=tpatterns{cvFolds(1,foldI)}; tpatternsTrain=double(tpatternsTrain');
        tpatternsTest=tpatterns{cvFolds(2,foldI)}; tpatternsTest=double(tpatternsTest');
        % define class lables
        labels=ones(72,1);
        labels(inanimates)=-1; % 1 = animate, -1 = inanimate
        % train and test the classifier
        model=svmtrain(labels,tpatternsTrain,libSVMsettings);
        [labelsPredicted,accuracy,decVals]=svmpredict(labels,tpatternsTest,model);
        accuracy_fold(foldI)=accuracy(1);
        
        % create null distribution for statistical inference
        for randomisationI=1:nRandomisations
            % randomise labels (for training)
            labelsRand=labels(randperm(length(labels)));
            % train and test the classifier using the randomised training labels
            modelRand=svmtrain(labelsRand,tpatternsTrain,libSVMsettings);
            [labelsPredictedRand,accuracyRand,decValsRand]=svmpredict(labels,tpatternsTest,modelRand);
            accuracy_randomisation_fold(randomisationI,foldI)=accuracyRand(1);
        end % randomisationI
    end % foldI
    
    % statistical inference
    accuracy=mean(accuracy_fold);
    accuracyH0=mean(accuracy_randomisation_fold,2);
    p=1-relRankIn_includeValue_lowerBound(accuracyH0,accuracy);
    
    % visualise results
    figure(10); clf;
    % plot null distribution
    hist(accuracyH0); hold on;
    % plot accuracy (mean across folds) found in the data
    xlim([5 95]); xlims=xlim;
    plot(accuracy,0,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',8);
    ylims=ylim;
    text(accuracy,0.04*ylims(2),'\bfdata','Color','r');
    % plot statistical result
    text(0.85*xlims(2),0.9*ylims(2),['p = ',sprintf('%1.4f',p)]);
    % label axes
    xlabel('accuracy');
    ylabel('frequency');
    title({'\fontsize{11}null distribution of classification accuracy',['\fontsize{8}',num2str(nRandomisations),' stimulus-label randomisations']})
end


%% RSA

if analyse.RSA
    % load behavioural RDM and image icons
    extra=load('pair1_subj1_extra.mat','behav_and_icons');
    
    judgmentRDM.RDM  = extra.behav_and_icons.RDM;
    judgmentRDM.name = 'similarity judgments';
    
    imageIcons  = extra.behav_and_icons.imageData;
    
    clear extra;
    
    % now we make an RDM per session
    RDMs=struct();
    for sessionI=1:2
        
        thisSess = sessions{sessionI};
        
        % the correlation distance patterns are computed using the pdist function
        RDMs(1,sessionI).RDM   = squareform(pdist(tpatterns{sessionI}','correlation'));
        RDMs(1,sessionI).name  = sprintf('hIT RDM | %s | session %d',thisSess,sessionI);
        RDMs(1,sessionI).color = [];
    end
    
    % show the 2 session RDMs
    figI=1;
    figure(figI);clf;set(gcf,'Position',[100 100 800 800],'Color','w')
    showRDMs(RDMs,figI);
    
    % ---------------------------------------------------------------------
    % compare RDMs across sessions 
    r12=corr([vectorizeRDMs(RDMs(1).RDM)]',[vectorizeRDMs(RDMs(2).RDM)]');
    % ---------------------------------------------------------------------
    
    avgRDM = averageRDMs_subjectSession(RDMs,'subject');
    avgRDM.name='hIT RDM averaged across sessions';
    
    figI=2;
    figure(figI);set(gcf,'Position',[100 100 800 800],'Color','w')
    showRDMs(avgRDM,figI);
    
    
    % define the labels and indices for familiar and unfamiliar images
    reductionLabels = {'familiar','unfamiliar'};
    reductionvectors = {1:36;37:72};
    nobjects = length(reductionvectors{1});
    
    reductionI=2;
    reduction = reductionvectors{reductionI};
    
    OwnBodyParts = 1:3; OwnFaces = 4:8; OwnPet = 9; OwnPlaces = 10:15; OwnObjects=16:18;
    OtherBodyParts = 19:21; OtherFaces = 22:26; OtherPet = 27; OtherPlaces = 28:33; OtherObjects=34:36;
    GeneralBodyParts = 37:44; GeneralFaces = 45:52; GeneralPets = [53 54]; GeneralPlaces=55:66; GeneralObjects=67:72;
    if reductionI==1
        bodies = [OwnBodyParts OtherBodyParts];%GeneralBodyParts-36;
        faces = [OwnFaces OwnPet OtherFaces OtherPet];%GeneralFaces GeneralPets]-36;
        places = [OwnPlaces OtherPlaces ];%GeneralPlaces-36;
        objects = [OwnObjects OtherObjects];%GeneralObjects-36;
    else
        bodies = GeneralBodyParts-36;
        faces =[GeneralFaces GeneralPets]-36;
        places = GeneralPlaces-36;
        objects = GeneralObjects-36;
    end
    nCols=4;
    cmap=RDMcolormap;
    colors=cmap([1 65 193 222],:);
    options.categoryColors=zeros(36,3);
    options.categoryColors(bodies,:)=repmat(colors(1,:),length(bodies),1);
    options.categoryColors(faces,:)=repmat(colors(2,:),length(faces),1);
    options.categoryColors(places,:)=repmat(colors(3,:),length(places),1);
    options.categoryColors(objects,:)=repmat(colors(4,:),length(objects),1);
    
    options.spheres=2;
    options.cols=options.categoryColors;
    options.replicability=0;
    options.view=1;
    
    D=avgRDM.RDM(reduction,reduction);
    [pats_mds_2D,stress,disparities]=extractMDS(D,2,options);
    
    % draw the mds
    nImages=size(pats_mds_2D,1);
    
    % compute image size
    imageAreaProportion=.5;
    boundingBoxArea=max(prod(range(pats_mds_2D)),max(range(pats_mds_2D))^2/10);
    totalImageArea=boundingBoxArea*imageAreaProportion;
    imageWidth=sqrt(totalImageArea/nImages);
    
    % smooth alpha channel
    transparentCol=[128 128 128 2];
    hsize=5*transparentCol(4);
    sigma=1*transparentCol(4);
    kernel=fspecial('gaussian', hsize, sigma);
    
    markerSize=85;
    figure(figI);clf;
    set(gcf,'Position',[ 100 100 800 800],'Color',[1 1 1],'Renderer','OpenGL','BackingStore','on'); % much better
    axes('Position',[0.05 0.2 0.9 0.75])
    hold on
    for imageI=1:nImages
        %[xs,ys,rgb3]=size(imageStruct(imageI).image);
        
        if reductionI==1
            transparent=imageIcons(imageI).image(:,:,1)==transparentCol(1) & imageIcons(imageI).image(:,:,2)==transparentCol(2) & imageIcons(imageI).image(:,:,3)==transparentCol(3);
        else
            transparent=imageIcons(imageI+36).image(:,:,1)==transparentCol(1) & imageIcons(imageI+36).image(:,:,2)==transparentCol(2) & imageIcons(imageI+36).image(:,:,3)==transparentCol(3);
        end
        if numel(transparentCol)==4
            % smooth alpha channel
            opacity=imfilter(double(1-transparent),kernel);
        else
            opacity=~transparent;
        end
        if reductionI==1
            if imageI<=18
                plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),...
                    'o','MarkerFaceColor',[128 128 128]./255,'MarkerEdgeColor',[128 128 128]./255,'MarkerSize',markerSize+20);
            end
        end
        
        %     plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),...
        %         'o','MarkerFaceColor',options.categoryColors(imageI,:),'MarkerEdgeColor',options.categoryColors(imageI,:),'MarkerSize',markerSize);
        %
        %imagesc(npats_mds_2D(imageI,1),npats_mds_2D(imageI,2),imageIcons(imageI+36).image);
        if reductionI==1
            image('CData',imageIcons(imageI).image,'XData',[pats_mds_2D(imageI,1)-imageWidth/2, pats_mds_2D(imageI,1)+imageWidth/2],'YData',[pats_mds_2D(imageI,2)+imageWidth/2, pats_mds_2D(imageI,2)-imageWidth/2],'AlphaData',opacity);
        else
            image('CData',imageIcons(imageI+36).image,'XData',[pats_mds_2D(imageI,1)-imageWidth/2, pats_mds_2D(imageI,1)+imageWidth/2],'YData',[pats_mds_2D(imageI,2)+imageWidth/2, pats_mds_2D(imageI,2)-imageWidth/2],'AlphaData',opacity);
        end
    end
    axis tight equal off
    annotation('textbox',[0 .90 1 0.1],'EdgeColor','none','String','MDS plot for unfamiliar images',...
        'HorizontalAlignment','center','FontWeight','bold','FontSize',18);
    
    axes('Position',[0.7 0 0.25 0.25])
    hold on;
    % add a micro mds plot
    markerSize=18;
    for imageI=1:nImages
        
        plot(pats_mds_2D(imageI,1),pats_mds_2D(imageI,2),...
            'o','MarkerFaceColor',options.categoryColors(imageI,:),'MarkerEdgeColor',options.categoryColors(imageI,:),'MarkerSize',markerSize);
    end
    axis tight equal off
      
    % relate hIT and judgments
    userOptions.analysisName = 'animacyVsJudgments';
    userOptions.candRDMdifferencesTest='conditionRFXbootstrap';
    candidateRDMs=cell(1); 
    candidateRDMs{1}=judgmentRDM;
    
    % ---------------------------------------------------------------------
    % compare the judgements RDM to an animacy model
    animates=[OwnBodyParts OwnFaces OwnPet OtherBodyParts OtherFaces OtherPet GeneralBodyParts GeneralFaces GeneralPets];
    inanimates=[OwnPlaces OwnObjects OtherPlaces OtherObjects GeneralPlaces GeneralObjects];
    
    model.RDM = ones(nImages*2,nImages*2);
    model.RDM(animates,animates) = 0;
    model.RDM(inanimates,inanimates) = 0;
    model.name = 'animacy';    
    
    candidateRDMs{2}= model;    
       
    stats_p_r=compareRefRDM2candRDMs(avgRDM.RDM, candidateRDMs, userOptions);
    
    % ---------------------------------------------------------------------
    
end

