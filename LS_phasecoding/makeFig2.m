% this script attempts to replicate figure 2.  It assumes the following:
%
% 1) you have buzcode (https://github.com/buzsakilab/buzcode) downloaded and in your matlab path
% 2) you are running this function from the dataset folder (with recording session subdirectories) or have changed the 'homeDirectory' variable 
%
%

%% plotting example cells
figure(2)

homeDirectory = pwd;

examples = {'20170509_468um_36um_170509_103451',                 2, 2,'ls';
%             '20170504_396um_0um_merge',                          19,1,'ls';
            '20170525_900um_936um_170525_150124' ,               8, 2,'ls';
            '20170509_468um_36um_170509_103451',                 24,1,'ls';
            '20170501_216um_0um_170501_112210',                  44,7,'ls';
            'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,1,'ls';
              };
          
for ex = 1:size(examples,1)
    cd(examples{ex,1})
    % load necessary data
    load([examples{ex,1} '.behavior.mat'])
    load([examples{ex,1} '.firingMaps.cellinfo.mat'])
    load([examples{ex,1} '.phaseMaps.cellinfo.mat'])
    load([examples{ex,1} '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
    
    spikes = bz_GetSpikes; % requires buzcode
    cond = examples{ex,3};
    cell = examples{ex,2};
    
    % start plotting
    subplot(5,5,1+((ex-1)*5))
    bz_plotTrials(behavior,'condition',cond,'spacing',6)
    % toggle the below line to show spike positions
%     scatter(phaseMaps.phaseMaps{cond}{cell}(:,3),phaseMaps.phaseMaps{cond}{cell}(:,4),'.m')
    axis square
    ylabel(['position (' behavior.units ')'])
    ylabel(['position (' behavior.units ')'])
    
    subplot(5,5,2+((ex-1)*5))
    imagesc(squeeze(firingMaps.rateMaps{cond}(cell,:,:)))
    ylabel('trial #')
    xlabel('linearized position')
    
    subplot(5,5,3+((ex-1)*5))
    if ~isempty(phaseMaps.phaseMaps{cond}{cell}) % prevents [] empty error
    scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),(phaseMaps.phaseMaps{cond}{cell}(:,end)),'.k');
    hold on
    scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),(phaseMaps.phaseMaps{cond}{cell}(:,end))+2*pi,'.k');
    nBins = length(behavior.events.map{1}.x);
    axis([0 nBins -pi pi*3])
    xlabel('linearized position')
    ylabel('theta phase')
    end
    
    subplot(5,5,4+((ex-1)*5))
    maxTau = max(positionDecodingMaxCorr_binned_box_mean.results{cell}.tau);
    rows = find(positionDecodingMaxCorr_binned_box_mean.results{cell}.condition==cond &...
        positionDecodingMaxCorr_binned_box_mean.results{cell}.tau==20);
    r = positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_rate(rows);
    p(:,1) = positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase_cos(rows);
    p(:,2) = positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase(rows);
    p(:,3) = positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase_all(rows);
    [t loc]=min(mean(p));
    p = p(:,loc);
    c = positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_chance_rate(rows);
    plot(ones(length(r),1)*2,r,'.r')
    hold on
    plot(ones(length(r),1)*3,p,'.g')
    plot(ones(length(r),1),c,'.k')
    errorbar(3,mean(p),std(p)*3,'g')
    errorbar(2,mean(r),std(r)*3,'r')
    errorbar(1,mean(c),std(c)*3,'k')
    axis([0 4 0 max([r' p' c'])*1.1])
    ylabel('mean squared error')
    xlabel('chance             rate               phase')
    [a rc] = ttest2(r,c);
    [a rp] = ttest2(r,p);
    [a pc] = ttest2(p,c);
    sigstar({[2,3],[1,2], [1,3]},[rp rc pc])
    title(['recording: ' examples{ex,1}])

%     saveFigure(['/home/david/Dropbox/Documents/pubs/ls_precession/' examples{ex,1} '_' num2str(cell)])
    
cd(homeDirectory)

clear p
end


%% now compile summary data...
hpc_phase = []; hpc_rate = [];
ls_rate = []; ls_phase = [];
ls_rec =[]; hpc_rec = [];
hpc_an = []; ls_an=[];
ls_depth=[];
hpc_depth=[];
ls_phase_info = []; hpc_phase_info =[];
ls_rate_info =[]; hpc_rate_info = [];
hpc_cell=[];
ls_cell=[];
ls_pvals = [];
hpc_shank = []; ls_shank = [];
hpc_field = []; ls_field = [];
hpc_reg = [];
ls_chance_rate = []; hpc_chance_rate = [];
ls_chance_phase = []; hpc_chance_phase = [];           

%%
hpc_pkFR = [];
hpc_stdev = [];
ls_pkFR = [];
ls_stdev = [];
ls_isi = [];
ls_waveform = [];
ls_meanRate = [];
%%
recordingList  = dir('*201*');

for i=1:length(recordingList)
   cd(recordingList(i).name) 
   if ~isempty(dir('*positionDecodingMaxCorr_binned_box_mean*')) & exist([recordingList(i).name '.placeFields.20_pctThresh.mat'])
        % load data for this session
        sessionInfo = bz_getSessionInfo;
        spikes = bz_GetSpikes;
        load([sessionInfo.FileName '.behavior.mat'])
        load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
        load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
        load([recordingList(i).name '.placeFields.20_pctThresh.mat'],'fields') 
        load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
        nBins = round(length(behavior.events.map{1}.x));
        conditions = length(unique(behavior.events.trialConditions));
   
        for cell =1:length(positionDecodingMaxCorr_binned_box_mean.results)
            if ~isempty(positionDecodingMaxCorr_binned_box_mean.results{cell})
            t_rate = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_rate',...
                'GroupingVariables',{'tau','condition'});
            t_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase',...
                'GroupingVariables',{'tau','condition'});
            t_phase_cos = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase_cos',...
                'GroupingVariables',{'tau','condition'});
            t_phase_chance = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_phase',...
                'GroupingVariables',{'tau','condition'});
            tab = join(join(join(t_rate,t_phase),t_phase_cos),t_phase_chance);
            for cond = 1:conditions
                if sum(behavior.events.trialConditions==cond) >= 10 % minimum number of trials
                if sum(sum(firingMaps.countMaps{cond}(cell,:,:))) >= 1.5 * sum(behavior.events.trialConditions==cond) % minimum firing rate
                %% carry on..
                rows = find(tab.condition==cond);
                kernelIDX = find(ismember(tab.tau(rows),20));
                [min_mse_rate] = tab.mean_mse_rate(rows(kernelIDX));
                [min_mse_phase_all] = tab.mean_mse_phase_cos(rows(kernelIDX));
                
                row = find(positionDecodingMaxCorr_binned_box_mean.results{cell}.condition == cond);
                col = find(positionDecodingMaxCorr_binned_box_mean.results{cell}.tau == 20);
                rows = intersect(row,col);
                [sig pval_vs_rate] = ttest2(positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase_cos(rows),positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_rate(rows));
                [sig_c pval_vs_control] = ttest2(positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase_cos(rows),positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_chance_phase(rows));
                additionalDepth = find(sessionInfo.spikeGroups.groups{spikes.shankID(cell)}==spikes.maxWaveformCh(cell))*10;
               if strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'hpc') | strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'ca3')  | strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'ca1') 
                   hpc_phase=[hpc_phase;min_mse_phase_all'];
                   hpc_rate=[hpc_rate;min_mse_rate'];
                   hpc_rec = [hpc_rec; i];
                   hpc_an = [hpc_an; sum(double(sessionInfo.animal))];
                   if ~isempty(sessionInfo.ca3)
                       hpc_reg = [hpc_reg; 3];
                   else
                       hpc_reg = [hpc_reg; 1];
                   end
                   if ~isempty(fields{cond}{cell}) 
                       hpc_field = [hpc_field;fields{cond}{cell}{1}.COM];
                   else
                       hpc_field = [hpc_field;nan];
                   end
                   if ~isempty(fields{cond}{cell})
                   hpc_pkFR = [hpc_pkFR; fields{cond}{cell}{1}.peakFR];
                   peakLoc = fields{cond}{cell}{1}.peakLoc;
                   hpc_stdev = [hpc_stdev;mean(std(firingMaps.rateMaps{cond}(cell,:,peakLoc),[],2))];
                   else
                   hpc_pkFR = [hpc_pkFR; nan];
                   hpc_stdev = [hpc_stdev;nan];
                   end
               elseif strcmp(positionDecodingMaxCorr_binned_box_mean.region{cell},'ls')    
                   if ~isempty(fields{cond}{cell})
                   ls_pkFR = [ls_pkFR; fields{cond}{cell}{1}.peakFR];
                   peakLoc = fields{cond}{cell}{1}.peakLoc;
                   ls_stdev = [ls_stdev;mean(std(firingMaps.rateMaps{cond}(cell,:,peakLoc),[],2))];
                   else
                   ls_pkFR = [ls_pkFR; nan];
                   ls_stdev = [ls_stdev;nan];
                   end
                   ls_waveform = [ls_waveform;spikes.rawWaveform{cell}];
                   ls_isi = [ls_isi; hist(diff(spikes.times{cell}),0:.001:.2)];
                   ls_meanRate = [ls_meanRate;sum(spikes.times{cell})./spikes.times{cell}(end)];
                   
                   ls_phase=[ls_phase;min_mse_phase_all'];
                   ls_rate=[ls_rate;min_mse_rate'];
                   ls_rec = [ls_rec; i];
                   ls_cell = [ls_cell; cell];
                   ls_an = [ls_an; sum(double(sessionInfo.animal))];
                   ls_pvals = [ls_pvals; pval_vs_control pval_vs_rate];
                   if ~isempty(fields{cond}{cell}) 
                       ls_field = [ls_field;fields{cond}{cell}{1}.COM];
                   else
                       ls_field = [ls_field;nan];
                   end
               end
               end
               end
            end
            end
        end
    end
cd(homeDirectory)
end

%% now plot summary data..

subplot(5,5,5)
histogram(ls_phase,2000:250:14000,'FaceColor','b')
hold on
histogram(ls_rate,2000:250:14000,'FaceColor','r')

subplot(5,5,10)
histogram(hpc_phase(~isnan(hpc_field))-hpc_rate(~isnan(hpc_field)),-6000:200:6000,'FaceColor','k','Normalization','pdf')
hold on
histogram(ls_phase-ls_rate,-6000:200:6000,'FaceColor','m','Normalization','pdf')

subplot(5,5,15)
for i=1:length(recordingList)
    f = find(hpc_rec==i);
    ff = find(~isnan(hpc_field));f = intersect(f,ff);
    hp(i) = nanmean(nanmean(hpc_phase(f,end)));
    hr(i) = nanmean(nanmean(hpc_rate(f,end)));
    f = find(ls_rec==i);
    if ~isempty(f)
        lr(i) = (nanmean(nanmean(ls_rate(f,end))));
        lp(i) = (nanmean(nanmean(ls_phase(f,end))));
    else
        lr(i) = nan;
        lp(i) = nan;
    end
end
histogram(hp-hr,-4000:180:4000,'FaceColor','k')
hold on
histogram(lp-lr,-4000:180:4000,'FaceColor','m')















