function [pats_mds_2D,stress,disparities,description]=extractMDS(D,nDims,options)
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here

import rsa.*
import rsa.fig.*
import rsa.fmri.*
import rsa.rdm.*
import rsa.sim.*
import rsa.spm.*
import rsa.stat.*
import rsa.util.*

options = setIfUnset(options,'MDScriterion','metricstress');
options = setIfUnset(options,'RDMname','RDM');
try
    [pats_mds_2D, stress, disparities] = mdscale(D, nDims,'criterion',options.MDScriterion);
    description{1}=sprintf('MDS computed using %s for %s',options.MDScriterion,options.RDMname);
catch
    try
        [pats_mds_2D, stress, disparities] = mdscale(D, nDims,'criterion','metricstress');
        description{1}=['\fontsize{12}MDS(reverted to stress, ',options.MDScriterion,' failed): ',options.RDMname];
    catch
        try
            D2=D+0.2;
            D2(logical(eye(length(D))))=0;
            [pats_mds_2D, stress, disparities] = mdscale(D2, nDims,'criterion',options.MDScriterion);
            description{1}=['\fontsize{12}MDS(added .2 to dissims to avoid colocalization)',options.RDMname];
        catch
            disp('MDS failed')
            description{1}='metricstress MDS failed...';
        end
    end
end
