function [varargout] = ConstructSearchlightRDMs(fullBrainVols,subject,userOptions)
%
% fMRISearchlight is a function which takes some full brain volumes of data,
% some binary mask and perfoms a searchlight in the data within
% each mask resulting in an RDM (lower triangular part) for each valid voxel mask location.
%
%
%       fullBrainVols --- The unmasked beta (or t) images.
%               fullBrainVols.(subject) is a [nVoxel nCondition nSession]-sized
%               matrix. The order of the voxels is that given by reshape or (:).
%
%        binaryMasks_nS --- The native- (subject-) space masks.
%               binaryMasks_nS.(subject).(mask) is a [x y z]-sized binary matrix
%               (the same size as the native-space 3D beta images.
%
%        userOptions --- The options struct.
%                userOptions.analysisName
%                        A string which is prepended to the saved files.
%                userOptions.rootPath
%                        A string describing the root path where files will be
%                        saved (inside created directories).
%                userOptions.subjectNames
%                        A cell array containing strings identifying the subject
%                        names. Defaults to the fieldnames in fullBrainVols.
%                userOptions.maskNames
%                        A cell array containing strings identifying the mask
%                        names. Defaults to the fieldnames of the first subject
%                        of binaryMasks_nS.
%                userOptions.voxelSize
%                        A tripple consisting of the [x y z] dimensions of each
%                        voxel in mm.
%                userOptions.structuralsPath
%                        A string which contains the absolute path to the
%                        location of the structural images and the normalisation
%                        warp definition file. It can contain the following
%                        wildcards which would be replaced as indicated:
%                                [[subjectName]]
%                                        To be replaced with the name of each
%                                        subject where appropriate.
%
%__________________________________________________________________________
% Ian Charest - October 2014
% Copyright (C) 2010 Medical Research Council


returnHere = pwd; % We'll come back here later

%% Set defaults and check options struct
if ~isfield(userOptions, 'analysisName'), error('fMRISearchlight:NoAnalysisName', 'analysisName must be set. See help'); end%if
if ~isfield(userOptions, 'rootPath'), error('fMRISearchlight:NoRootPath', 'rootPath must be set. See help'); end%if
userOptions = setIfUnset(userOptions, 'subjectNames', fieldnames(fullBrainVols));
if ~isfield(userOptions, 'voxelSize'), error('fMRISearchlight:NoVoxelSize', 'voxelSize must be set. See help'); end%if

% The analysisName will be used to label the files which are eventually saved.
RDMsFilename = [userOptions.analysisName,'_',subject,'_Searchlight_RDMs.mat'];
DetailsFilename = [userOptions.analysisName,'_',subject,'_Searchlight_Details.mat'];

promptOptions.functionCaller = 'ConstructSearchlightRDMs';
promptOptions.defaultResponse = 'S';
promptOptions.checkFiles(1).address = fullfile(userOptions.rootPath, 'RDMs', RDMsFilename);
promptOptions.checkFiles(2).address = fullfile(userOptions.rootPath, 'Details', DetailsFilename);

overwriteFlag = overwritePrompt(userOptions, promptOptions);

voxSize_mm = userOptions.voxelSize;
searchlightRad_mm = userOptions.searchlightRadius;	

if overwriteFlag
	
	searchlightOptions.monitor = false;
	searchlightOptions.fisher = true;
	
	fprintf('***\t Shining RSA searchlights...\t***\n');
    fprintf('***\t Constructing searchlights RDMs for subject %s ...\t***\n',subject);

    
    t_pats = fullBrainVols.(subject)'; clear fullBrainVols;
    
    nConditions = size(t_pats, 1);
    nComps = ((nConditions*nConditions)-nConditions)/2;
    
    % Prepare masks
    readPath = replaceWildcards(userOptions.maskPath, '[[subjectName]]', subject,'[[maskName]]', userOptions.maskNames{1});
    mask = spm_read_vols(spm_vol(readPath));
            
    
    mask(isnan(mask)) = 0; % Just in case!
    if ndims(mask)==3
        inputDataMask=logical(mask);
        mappingMask_request=logical(mask);
    else
        inputDataMask=logical(mask(:,:,:,1));
        mappingMask_request=logical(mask(:,:,:,2));
    end
        
    if (size(t_pats,2)>sum(inputDataMask(:)))
			t_pats=t_pats(:,inputDataMask(:));
    end%if
    
	% Other data
	volSize_vox=size(inputDataMask);
	rad_vox=searchlightRad_mm./voxSize_mm;
	minMargin_vox=floor(rad_vox);


	%% create spherical multivariate searchlight
	[x,y,z]=meshgrid(-minMargin_vox(1):minMargin_vox(1),-minMargin_vox(2):minMargin_vox(2),-minMargin_vox(3):minMargin_vox(3));
	sphere=((x*voxSize_mm(1)).^2+(y*voxSize_mm(2)).^2+(z*voxSize_mm(3)).^2)<=(searchlightRad_mm^2);  % volume with sphere voxels marked 1 and the outside 0
	sphereSize_vox=[size(sphere),ones(1,3-ndims(sphere))]; % enforce 3D (matlab stupidly autosqueezes trailing singleton dimensions to 2D, try: ndims(ones(1,1,1)). )

	% compute center-relative sphere SUBindices
	[sphereSUBx,sphereSUBy,sphereSUBz]=ind2sub(sphereSize_vox,find(sphere)); % (SUB)indices pointing to sphere voxels
	sphereSUBs=[sphereSUBx,sphereSUBy,sphereSUBz];
	ctrSUB=sphereSize_vox/2+[.5 .5 .5]; % (c)en(t)e(r) position (sphere necessarily has odd number of voxels in each dimension)
	ctrRelSphereSUBs=sphereSUBs-ones(size(sphereSUBs,1),1)*ctrSUB; % (c)en(t)e(r)-relative sphere-voxel (SUB)indices

	
	%% define masks
	validInputDataMask=inputDataMask;

	sumAbsY=sum(abs(t_pats),1);

	validYspace_logical= (sumAbsY~=0) & ~isnan(sumAbsY); clear sumAbsY;
	validInputDataMask(inputDataMask)=validYspace_logical; % define valid-input-data brain mask

	t_pats=t_pats(:,validYspace_logical); % reduce t_pats to the valid-input-data brain mask
	nVox_validInputData=size(t_pats,2);
	
	mappingMask_request_INDs=find(mappingMask_request);
	nVox_mappingMask_request=length(mappingMask_request_INDs);

	volIND2YspaceIND=nan(volSize_vox);
	volIND2YspaceIND(validInputDataMask)=1:nVox_validInputData;

	% n voxels contributing to infobased t at each location
	n=nan(1,nVox_mappingMask_request);

	%% similarity-graph-map the volume with the searchlight
	
	searchlightRDMs = nan(nVox_mappingMask_request,nComps);
	
	%% THE BIG LOOP! %%

	for cMappingVoxI=1:nVox_mappingMask_request
		
		if mod(cMappingVoxI,1000)==0			
			fprintf('.');			
		end%if

		[x,y,z]=ind2sub(volSize_vox,mappingMask_request_INDs(cMappingVoxI));

		% compute (sub)indices of (vox)els (c)urrently (ill)uminated by the spherical searchlight
		cIllVoxSUBs=repmat([x,y,z],[size(ctrRelSphereSUBs,1) 1])+ctrRelSphereSUBs;

		% exclude out-of-volume voxels
		outOfVolIs=(cIllVoxSUBs(:,1)<1 | cIllVoxSUBs(:,1)>volSize_vox(1)|...
					cIllVoxSUBs(:,2)<1 | cIllVoxSUBs(:,2)>volSize_vox(2)|...
					cIllVoxSUBs(:,3)<1 | cIllVoxSUBs(:,3)>volSize_vox(3));

		cIllVoxSUBs=cIllVoxSUBs(~outOfVolIs,:);

		% list of (IND)ices pointing to (vox)els (c)urrently (ill)uminated by the spherical searchlight
		cIllVox_volINDs=sub2ind(volSize_vox,cIllVoxSUBs(:,1),cIllVoxSUBs(:,2),cIllVoxSUBs(:,3));

		% restrict searchlight to voxels inside validDataBrainMask
		cIllValidVox_volINDs=cIllVox_volINDs(validInputDataMask(cIllVox_volINDs));
		cIllValidVox_YspaceINDs=volIND2YspaceIND(cIllValidVox_volINDs);

		% note how many voxels contributed to this locally multivariate stat
		n(cMappingVoxI)=length(cIllValidVox_YspaceINDs);
		
		if n(cMappingVoxI) < 2, continue; end%if % This stops the function crashing if it accidentally encounters an out-of-brain floating voxel (these can occur if, for example, skull stripping fails)
		
		% Locally store the full brain's worth of indexed RDMs.
		searchlightRDMs(cMappingVoxI,:) = pdist(t_pats(:,cIllValidVox_YspaceINDs), 'correlation');
		
	end%for:cMappingVoxI
    
    %% Save relevant info

	timeStamp = datestr(now);

	fprintf(['Saving searchlight RDMs to ' fullfile(userOptions.rootPath, 'RDMs', RDMsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'RDMs');
	save(RDMsFilename, 'searchlightRDMs','mappingMask_request_INDs','volSize_vox','-v7.3');
	
	fprintf(['Saving Details to ' fullfile(userOptions.rootPath, 'Details', DetailsFilename) '\n']);
	gotoDir(userOptions.rootPath, 'Details');
	save(DetailsFilename, 'timeStamp', 'userOptions');
    
else
    load(fullfile(userOptions.rootPath, 'RDMs',RDMsFilename), 'searchlightRDMs','mappingMask_request_INDs','volSize_vox');
    
end%if


if nargout == 3
	varargout{1} = searchlightRDMs;
	varargout{2} = mappingMask_request_INDs;
	varargout{3} = volSize_vox;
elseif nargout > 0
	error('0, or 3 arguments out, please.');
end%if:nargout

cd(returnHere); % And go back to where you started

end%function

