


recordingList  = dir('*201*');
homeDirectory = pwd;
figure(4)
clf
ls_rate = []; ls_phase = [];
ls_tau_phase = []; ls_tau_rate = [];
ls_phase_pval = [];
ls_rate_pval = [];
ls_rec =[]; 
ls_an=[];
ls_depth=[];
ls_cell=[];
ls_shank = [];

for i=1:length(recordingList)
    cd(recordingList(i).name) 
    if ~isempty(dir('*positionDecodingMaxCorr_binned_box_mean*')) & exist([recordingList(i).name '.placeFields.20_pctThresh.mat'])
        sessionInfo = bz_getSessionInfo(pwd,'noprompts',true);
        spikes = bz_GetSpikes;
        load([sessionInfo.FileName '.behavior.mat'])
        load([recordingList(i).name '.firingMaps.cellinfo.mat'],'firingMaps')
        load([recordingList(i).name '.phaseMaps.cellinfo.mat'],'phaseMaps')
        [binnedfiringMaps.phaseMaps] = bz_phaseMap2Bins(phaseMaps.phaseMaps,firingMaps.rateMaps,behavior);
        load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
        conditions = length(unique(behavior.events.trialConditions));
            
        for cell =1:length(positionDecodingMaxCorr_binned_box_mean.results)
            if ~isempty(positionDecodingMaxCorr_binned_box_mean.results{cell})
            t_rate = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_rate',...
                'GroupingVariables',{'tau','condition'});
            t_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase',...
                'GroupingVariables',{'tau','condition'});
            t_phase_cos = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase_cos',...
                'GroupingVariables',{'tau','condition'});
            t_phase_all = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase_all',...
                'GroupingVariables',{'tau','condition'});
            t_chance_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_rate',...
                'GroupingVariables',{'tau','condition'});
            t_chance_rate = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_phase',...
                'GroupingVariables',{'tau','condition'});
            tab = join(join(join(join(join(t_rate,t_phase),t_chance_rate),t_chance_phase),t_phase_cos),t_phase_all);
            for cond = 1:conditions
               if sum(behavior.events.trialConditions==cond) >= 10 %%%%%%%%%%%%%%%%%%%%%%%%%%
               if sum(sum(firingMaps.countMaps{cond}(cell,:,:))) >= 1.5 * sum(behavior.events.trialConditions==cond) 
                    nTrials = sum(behavior.events.trialConditions==cond);
                %% carry on..
                rows = find(tab.condition==cond);    
                kernelIDX = find(ismember(tab.tau(rows),20));
                [min_mse_rate] = tab.mean_mse_rate(rows(kernelIDX));
                [min_mse_phase_all] = tab.mean_mse_phase_cos(rows(kernelIDX));
               if strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'hpc') | strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'ca3')  | strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'ca1') 
               % do nothing
               elseif strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'ls') 
                   ls_phase=[ls_phase;min_mse_phase_all'];
                   ls_rate=[ls_rate;min_mse_rate'];
                   ls_rec = [ls_rec; i];
                   ls_an = [ls_an; sessionInfo.animal];
                   additionalDepth = find(sessionInfo.spikeGroups.groups{spikes.shankID(cell)}==spikes.maxWaveformCh(cell))*10;
                   ls_depth = [ls_depth; (sessionInfo.depth)+additionalDepth];
                   ls_cell = [ls_cell; cell];
                   ls_shank = [ls_shank;spikes.shankID(cell)];
               end
               end
               end
            end
        end
        end
    end
    cd(homeDirectory)
end

for i=1:length(recordingList)
f = find(ls_rec==i);
if ~isempty(f)
anim(i) = sum(double(ls_an(f(1),:)));
depths(i) = ls_depth(f(1));
end
end
u = unique(anim);
for i=1:length(recordingList)
    f = find(ls_rec==i);
    if ~isempty(f)
        lr(i) = (nanmean(nanmean(ls_rate(f,end))));
        lp(i) = (nanmean(nanmean(ls_phase(f,end))));
    else
    end
end

offsets = [0 500 200 0 2500 600];
for i=1:length(u)
    f = find(anim==u(i));
    ff=find(sum(double(ls_an)')==u(i));
    if ~isempty(ff)
        if i == 2
        %     f = f(1:end-3);
        elseif i == 3
        %     f = f(1:end-3);
        elseif i == 5
            f = f(1:end-1);
        end
        subplot(3,2,i);hold on
        plot((lp(f))-(lr(f)),-depths(f),'.k')
        [a  b] = corr(depths(f)',(lp(f)')-(lr(f))');
        [x y] = polyfit(depths(f),(lp(f))-(lr(f)),1);
        y1 = polyval(x,depths(f));
        plot(y1,-depths(f),'r')
        title([ls_an(ff(1),:) ': r = ' num2str(a)])
        axis([-3000 3000 min(-depths(f)) max(-depths(f))])
    end
end



