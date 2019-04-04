clc
clear

% this is heavily inspired by the great MEG MVPA tutorial of
% Matthias Guggenmos
% https://github.com/m-guggenmos/megmvpa

load('/media/rdsfs2/irsaeeg_2016/Data_old/YC_180787/SET-ERPs/mat_binIMG_mrk_p_YC_180787.mat')

desiredtimes = 154:length(times);

workingdir = pwd;

addpath(genpath(fullfile(workingdir,'software','GLMdenoise-1.4')));
addpath(genpath(fullfile(workingdir,'software','rsatoolbox')));
addpath(genpath(fullfile(workingdir,'software','libsvm-3.22')));

% rearrange epochs, and select scalp electrodes
epochs = permute(catcell(3,C),[3 1 2]);
epochs = epochs(:,1:128,desiredtimes);

times = times(desiredtimes);

labels = repmat(1:72,30,1); labels = labels(:);
triggers = unique(labels);

% Parameters
n_perm = 20;  % number of permutations
n_pseudo = 10;  % number of pseudo-trials
n_conditions = length(triggers);
n_sensors = size(epochs, 2);
n_time = size(epochs, 3);

% n_sessions = length(sessions);
% preallocate results
result_cv = nan(n_perm, n_conditions, n_conditions, n_time);

X = epochs;
y = labels;

conditions = unique(y);
n_trials = histc(y, conditions);

%% now we compute the RDMs

for f = 1:n_perm
    fprintf('\tPermutation %g / %g\n', f, n_perm)
    
    % precompute permutations
    ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    ind_pseudo_test = nan(n_conditions, n_conditions, 2);
    labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
    labels_pseudo_test = nan(n_conditions, n_conditions, 2);
    for c1 = 1:n_conditions
        range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for c2 = 1:n_conditions
            range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
            ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
            ind_pseudo_test(c1, c2, :) = [c1 c2];
            labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                [conditions(c1)*ones(1,n_pseudo-1) conditions(c2)*ones(1,n_pseudo-1)];
            labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
        end
    end
    train_indices = cell(1, n_conditions*(n_pseudo-1));
    test_indices = cell(1, n_conditions);
    for c1 = 1:n_conditions  % separate permutation for each condition
        prm_ = randperm(n_trials(c1));
        prm = cell(1, n_pseudo);
        splitsize = n_trials(c1) / n_pseudo;
        for i = 1:n_pseudo
            idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
            prm{i} = prm_(idxs + 1);
        end
        ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
        xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
        for i = 1:length(xrange)
            train_indices{xrange(i)} = ind{i};
        end
        test_indices{c1} = ind{end};
    end
    
    % 1. Compute pseudo-trials for training and test
    Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
    Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
    for i = 1:length(train_indices)
        Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
    end
    for i = 1:length(test_indices)
        Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
    end
    
    
    % 2. Whitening using the Epoch method
    sigma_conditions = reshape(squeeze(labels_pseudo_train(1,:,n_pseudo:end))',1,[]);
    sigma_ = nan(n_conditions, n_sensors, n_sensors);
    for c = 1:n_conditions
        % compute sigma for each time point, then average across time
        tmp_ = nan(n_time, n_sensors, n_sensors);
        for t = 1:n_time
            tmp_(t, :, :) = rsa.stat.cov1para(Xpseudo_train(sigma_conditions==c, :, t));
        end
        sigma_(c, :, :) = mean(tmp_, 1);
    end
    sigma = squeeze(mean(sigma_, 1));  % average across conditions
    sigma_inv = sigma^-0.5;
    for t = 1:n_time
        Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
        Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
    end
    
    for t = 1:n_time
        for c1 = 1:n_conditions-1
            for c2 = c1+1:n_conditions
                % 3. Fit the classifier using training data
                data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                y_train = squeeze(labels_pseudo_train(c1, c2, :));
                model_svm = svmtrain(y_train, data_train, '-c 1 -q 0 -t 0');
                
                % 4. Compute and store classification accuracies
                data_test = Xpseudo_test(squeeze(ind_pseudo_test(c1, c2, :)), :, t);
                y_test = squeeze(labels_pseudo_test(c1, c2, :));
                result_cv(f, c1, c2, t) = ...
                    mean(svmpredict(y_test,data_test,model_svm,'-q 0 -t 0')==y_test)-0.5;
                
            end
        end
    end
end
result_cv_ = squeeze(nanmean(result_cv, 1));
save(fullfile(workingdir, 'result_svm_cv.mat'), 'result_cv','result_cv_');

%% visualise
figure; plot(times,squeeze(nanmean(nanmean(result_cv_,1),2)));

% show the RDM as a movie
figure;
for timeI=1:20:length(times)
    thistime = times(timeI);
    thisrdm = result_cv_(:,:,timeI);
    imagesc(thisrdm);
    colormap(rsa.fig.RDMcolormap);
    colorbar;
    title(sprintf('time: %.2f',thistime))
    pause(.25)
end


