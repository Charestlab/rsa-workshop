function [conditionVec,partition] = getSPMconditionVec(SPM,conditionLabels);

% getSPMconditionVec is a function which will extract from the condition vector 
% from the SPM structure and define the partitioning.
%
% function [conditionVec,partition] = getSPMconditionVec(SPM,Opt.conditionLabels);
%
%        SPM --- The SPM structure.
%
%        conditionLabels --- The condition labels.
%                these are also stored in {SPM.Sess(1).Fc.name}   
%                
%  
%  Ian Charest 05-01-2017

import rsa.*
import rsa.spm.*
import rsa.util.*

nSessions   = length(SPM.Sess);
nVolumes    = numel(SPM.Sess(1).row);
nConditions = length(conditionLabels); 
% verify that all sessions have the same number of volumes
for sessionI = 1:nSessions
    assert(numel(SPM.Sess(sessionI).row)==nVolumes,sprintf('error: not all sessions have %d volumes - check nVolumes for session %d\n',nVolumes,sessionI));
end

% define the condition vector
nBetas =  numel(SPM.Vbeta);
conditionVec = zeros(nBetas,1);
for betaI = 1:nBetas
    for conditionI=1:nConditions
        if strfind(SPM.Vbeta(betaI).descrip,conditionLabels{conditionI});
            conditionVec(betaI) = conditionI;
        end            
    end    
end

% define the partition vector
partitionlist = 1:nSessions;
partition = zeros(1,nBetas);
for conditionI = 1:nConditions
    partition(conditionVec==conditionI) = partitionlist;
end







