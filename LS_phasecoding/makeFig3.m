

clf
l_rate_maps_smooth = [];
l_phase_maps_smooth = [];
h_phase_maps_smooth = [];
h_rate_maps_smooth = [];
hpc_field = [];
ls_field = [];
h_an = [];
l_an =[];
recordingList  = dir('*201*');
homeDirectory = pwd;
tau = 20; % smoothing window (# of bins)


%% get rate/phase maps from all sessions
for i=1:length(recordingList)
   cd(recordingList(i).name) 
    if ~isempty(dir('*positionDecodingMaxCorr_binned_box_mean*')) & exist([recordingList(i).name '.placeFields.20_pctThresh.mat'])
        sessionInfo = bz_getSessionInfo(pwd,'noprompts',true);
        animal = sessionInfo.animal;
        spikes = bz_GetSpikes;
        load([sessionInfo.FileName '.behavior.mat'])
         [firingMaps] = bz_firingMap1D(spikes,behavior,tau);
        load([recordingList(i).name '.phaseMaps.cellinfo.mat'],'phaseMaps')
        load([recordingList(i).name '.placeFields.20_pctThresh.mat'],'fields') 
        [binnedfiringMaps.phaseMaps] = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
        conditions = length(unique(behavior.events.trialConditions));
        for cell =1:length(spikes.times)
            for cond = 1:conditions
               if sum(behavior.events.trialConditions==cond) >= 10 %%%%%%%%%%%%%%%%%%%%%%%%%%
               if sum(sum(firingMaps.countMaps{cond}(cell,:,:))) >= 1.5 * sum(behavior.events.trialConditions==cond) 
               if strcmp(spikes.region{cell},'hpc') | strcmp(spikes.region{cell},'ca3')  | strcmp(spikes.region{cell},'ca1') 
                    for field = 1:length(fields{cond}{cell})
                        if ~isempty(fields{cond}{cell}) %& length(fields{cond}{cell}) == 1
                           hpc_field = [hpc_field;fields{cond}{cell}{field}.COM];
                        else
                           hpc_field = [hpc_field;nan];
                        end
                        hh = (squeeze((squeeze(binnedfiringMaps.phaseMaps{cond}(cell,:,:)))));
                        for trial = 1:size(hh,1)
                        hh(trial,:) = circ_smoothTS(hh(trial,:),tau,'method','mean','exclude',0); 
                        end
                        h_rate_maps_smooth = [h_rate_maps_smooth;(squeeze(mean(firingMaps.rateMaps_box{cond}(cell,:,:),2))')]; 
                        h_phase_maps_smooth = [h_phase_maps_smooth;(squeeze(circ_mean(hh)))]; clear hh
                        h_an = [h_an; sum(double(animal))];
                    end
                    if isempty(fields{cond}{cell})
                        hpc_field = [hpc_field;nan];
                        hh = (squeeze((squeeze(binnedfiringMaps.phaseMaps{cond}(cell,:,:)))));
                        for trial = 1:size(hh,1)
                        hh(trial,:) = circ_smoothTS(hh(trial,:),tau,'method','mean','exclude',0); 
                        end
                        h_rate_maps_smooth = [h_rate_maps_smooth;(squeeze(mean(firingMaps.rateMaps_box{cond}(cell,:,:),2))')]; 
                        h_phase_maps_smooth = [h_phase_maps_smooth;(squeeze(circ_mean(hh)))]; clear hh
                        h_an = [h_an; sum(double(animal))];
                    end
               elseif strcmp(spikes.region{cell},'ls') 
                   field = 1;%for field = 1:length(fields{cond}{cell})
                    if ~isempty(fields{cond}{cell}) %& length(fields{cond}{cell}) == 1
                        ls_field = [ls_field;fields{cond}{cell}{field}.COM];
                    else
                        ls_field = [ls_field;nan];
                    end
                    ll = (squeeze((squeeze(binnedfiringMaps.phaseMaps{cond}(cell,:,:)))));
                    for trial = 1:size(ll,1)
                       ll(trial,:) = circ_smoothTS(ll(trial,:),tau,'method','mean','exclude',0); 
                    end
                    l_rate_maps_smooth = [l_rate_maps_smooth;squeeze(mean(firingMaps.rateMaps_box{cond}(cell,:,:),2))'];
                    l_phase_maps_smooth = [l_phase_maps_smooth;(circ_mean(ll))];clear ll
                    l_an = [l_an; sum(double(animal))];
               end
               end
               end
            end
        end
    end
cd(homeDirectory)
end


%% plotting 
figure(3)

subplot(4,4,1)
imagesc(corr(h_rate_maps_smooth))
caxis([.1 1])
title('HPC rate (all)')

subplot(4,4,2)
imagesc(corr(h_rate_maps_smooth(~isnan(hpc_field),:)))
caxis([.1 1])
title('HPC rate (place fields)')

subplot(4,4,3)
imagesc(corr(l_rate_maps_smooth))
caxis([.1 1])
title('LS rate (all)')

subplot(4,4,4)
imagesc(corr(l_rate_maps_smooth(~isnan(ls_field),:)))
caxis([.1 1])
title('LS rate (place fields)')

subplot(4,4,5)
for i=1:size(h_phase_maps_smooth,2)
    for j = 1:size(h_phase_maps_smooth,2)
        [h_phase_corr_all(i,j) p] = circ_corrcc(h_phase_maps_smooth(:,i),h_phase_maps_smooth(:,j));
    end
end
imagesc(h_phase_corr_all)
caxis([.1 1])
title('HPC phase (all)')

subplot(4,4,6)
for i=1:size(h_phase_maps_smooth,2)
    for j = 1:size(h_phase_maps_smooth,2)
        [h_phase_corr(i,j) p] = circ_corrcc(h_phase_maps_smooth(~isnan(hpc_field),i),h_phase_maps_smooth(~isnan(hpc_field),j));
    end
end
imagesc(h_phase_corr)
caxis([.1 1])
title('HPC phase (place fields)')

subplot(4,4,7)
for i=1:size(l_phase_maps_smooth,2)
    for j = 1:size(l_phase_maps_smooth,2)
        [l_phase_corr(i,j) p] = circ_corrcc(l_phase_maps_smooth(:,i),l_phase_maps_smooth(:,j));
    end
end
imagesc(l_phase_corr)
caxis([.1 1])
title('LS phase (all)')

subplot(4,4,8)
for i=1:size(l_phase_maps_smooth,2)
    for j = 1:size(l_phase_maps_smooth,2)
        [l_phase_corr(i,j) p] = circ_corrcc(l_phase_maps_smooth(~isnan(ls_field),i),l_phase_maps_smooth(~isnan(ls_field),j));
    end
end
imagesc(l_phase_corr)
caxis([.1 1])
title('LS phase (place fields)')

subplot(4,4,9)
hold on
plot(diag(flipud(h_phase_corr_all)),'b')
plot(diag(flipud(corr(h_rate_maps_smooth))),'r')
axis([0 200 0 1])

subplot(4,4,10)
hold on
plot(diag(flipud(h_phase_corr)),'b')
plot(diag(flipud(corr(h_rate_maps_smooth(~isnan(hpc_field),:)))),'r')
axis([0 200 0 1])

subplot(4,4,11)
hold on
plot(diag(flipud(l_phase_corr)),'b')
plot(diag(flipud(corr(l_rate_maps_smooth))),'r')
axis([0 200 0 1])

subplot(4,4,12)
hold on
plot(diag(flipud(l_phase_corr)),'b')
plot(diag(flipud(corr(l_rate_maps_smooth(~isnan(ls_field),:)))),'r')
axis([0 200 0 1])






