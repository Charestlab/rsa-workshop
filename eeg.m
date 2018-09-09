dataroot = '/media/charesti-start/data/irsa-eeg/BIDS/derivatives/eegprep/';
% trials x channels x timepoints
data = load(fullfile(dataroot, 'sub-01', 'sub-01_task-irsa_epo.mat'));

% plot grand average per trial
plot(squeeze(mean(data.epochs, 1)).');

% conditions = unique(data.events);
% dims = size(data.epochs);
% nconditions = numel(conditions);
% nchannels = dims(2);
% ntimepoints = dims(3);
% data_reshaped = zeros([nchannels ntimepoints nconditions ntrials??]);
% %reshape to electrodes x time x conditions x trials
% for c=1:length(conditions)
%     disp(c);
% end
% %B = permute(A,[3 2 1]);


