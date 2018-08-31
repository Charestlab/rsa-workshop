function [pats_mds_2D,stress,disparities]=extractMDS(D,nDims,options)


try
    [pats_mds_2D, stress, disparities] = mdscale(D, nDims,'criterion',options.MDScriterion);
catch
    try
        [pats_mds_2D, stress, disparities] = mdscale(D, nDims,'criterion','metricstress');
        %description{1}=['\fontsize{12}MDS(reverted to stress, ',options.MDScriterion,' failed): ',RDMname];
    catch
        try
            D2=D+0.2;
            D2(logical(eye(length(D))))=0;
            [pats_mds_2D, stress, disparities] = mdscale(D2, nDims,'criterion',options.MDScriterion);
            %description{1}=['\fontsize{12}MDS(added .2 to dissims to avoid colocalization)',RDMname];
        catch
            disp('MDS failed')
        end
    end
end