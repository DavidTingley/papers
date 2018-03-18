% replicates figure 5 from the paper


homeDirectory = pwd;
%             recording                            hpc      ls  condition
examples = {'20161018_1552um_2304um_merge',           70,     10, 2;...
            '20170501_216um_0um_170501_112210',     106,    3, 4;...
            '20170505_396um_0um_merge',             63,     6, 4;...
            '20170506_468um_36um_170506_111835',    105,    18, 5;...
            '20170509_468um_36um_170509_103451',    72,     5, 5};
        
%% plot example cell pairs
for ex = 1:length(examples)
   
   cd(examples{ex,1})
   sessionInfo = bz_getSessionInfo;
   load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
   load([sessionInfo.FileName '.behavior.mat'])
   cond = examples{ex,4};
   hpc = examples{ex,2};
   ls = examples{ex,3};
   
   subplot(6,5,ex)
   imagesc(squeeze(firingMaps.rateMaps{cond}(hpc,:,:)))
   
   subplot(6,5,ex+5)
   scatter(phaseMaps.phaseMaps{cond}{hpc}(:,1),(phaseMaps.phaseMaps{cond}{hpc}(:,end)),'.k');
   hold on
   scatter(phaseMaps.phaseMaps{cond}{hpc}(:,1),(phaseMaps.phaseMaps{cond}{hpc}(:,end))+2*pi,'.k');
   axis([0 200 -pi pi*3])
   
   subplot(6,5,ex+10)
   scatter(phaseMaps.phaseMaps{cond}{ls}(:,1),(phaseMaps.phaseMaps{cond}{ls}(:,end)),'.k');
   hold on
   scatter(phaseMaps.phaseMaps{cond}{ls}(:,1),(phaseMaps.phaseMaps{cond}{ls}(:,end))+2*pi,'.k');
   axis([0 200 -pi pi*3])
   
   subplot(6,5,ex+15)
   spikes = bz_GetSpikes;
   for spk = 1:length(spikes.times) % restrict to only spike times during this block of trials
      spikes.times{spk} = Restrict(spikes.times{spk},behavior.events.trialIntervals(behavior.events.trialConditions==cond,:));
   end
   [times groups] = spikes2sorted(spikes.times);
   [ccg ts] = CCG(times,groups,'binSize',.001,'duration',1);
   plot(ts,smooth(ccg(:,ls,hpc),20)*20,'k')
   axis([-.5 .5 0 max(smooth(ccg(:,hpc,ls),20))*20])
   
   subplot(6,5,ex+20)
    try
    load('assembliesCrossRegion_split_w_theta_08-Nov-2017.mat','dev*','pairs','coords');
    catch
    load('assembliesCrossRegion_split_w_theta.mat','dev*','pairs');%c
    end
    f = find(pairs(:,2) == hpc); ff = find(pairs(:,1)==ls);
    pair = intersect(f,ff);
    plot(dev{cond}(:,pair),'k')
    hold on
    plot(nanmean(devControl{cond}(:,pair,:),3),'r')
    axis([0 200 min(dev{cond}(:,pair))-std(dev{cond}(:,pair)) max(dev{cond}(:,pair))+std(dev{cond}(:,pair))])
    [a b] =  min(dev{cond}(:,pair));
    [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
    imp = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));   
    title(['assembly strength:' num2str(imp)])
    cd(homeDirectory)
end

%% compile dataset
        % compiling cell assembly data across regions and recordings...
recordingList  = dir('*201*');
assemblyStrength = [];
optimalWindow = [];
phaseRateBias = [];
noAssembly_phaseRateBias = [];
pairCount = 0;

for i=1:length(recordingList)
    cd(recordingList(i).name) 
    if exist('assembliesCrossRegion_split_w_theta_08-Nov-2017.mat') || exist('assembliesCrossRegion_split_w_theta.mat')
%         try
%         load('assembliesCrossRegion_split_w_theta_08-Nov-2017.mat','dev*','pairs','coords');
%         catch
%         load('assembliesCrossRegion_split_w_theta.mat','dev*','pairs');
%         end
        try
        load('assembliesCrossRegion_split_w_theta.mat','dev*','pairs');
        catch
        load('assembliesCrossRegion_split_w_theta_08-Nov-2017.mat','dev*','pairs','coords');%c
        end
        sessionInfo =bz_getSessionInfo;
        load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
        conditions = length(unique(positionDecodingMaxCorr_binned_box_mean.results{1}.condition));
        for cell =1:length(positionDecodingMaxCorr_binned_box_mean.results)
            t_rate = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_rate',...
                'GroupingVariables',{'tau','condition'});
            t_phase_cos = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase_cos',...
                'GroupingVariables',{'tau','condition'});
            t_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase',...
                'GroupingVariables',{'tau','condition'});
            t_chance_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_phase',...
                'GroupingVariables',{'tau','condition'});
            tab = join(join(join(t_rate,t_phase),t_chance_phase),t_phase_cos);
            for cond = 1:conditions
                rows = find(tab.condition==cond);
                tabIndex = find(ismember(tab.tau(rows),20));
                if cond <= length(dev)  && sum(behavior.events.trialConditions==cond) > 10 && ~isempty(tabIndex)               
                    min_mse_rate(cell,cond) = min(tab.mean_mse_rate(rows(tabIndex)));
                    min_mse_phase_all(cell,cond) = min([tab.mean_mse_phase_cos(rows(tabIndex)),tab.mean_mse_phase(rows(tabIndex))]);
                else
                    min_mse_rate(cell,cond) = nan;
                    min_mse_phase_all(cell,cond) = nan;
                end
            end
        end
        
        sessionInfo = bz_getSessionInfo;
        spikes = bz_GetSpikes('noprompt',true);
        if ~isempty(pairs)
            for cond = 1:length(dev)
               p=[];
               pp=[];
               if sum(behavior.events.trialConditions==cond) > 10
               for pair = 1:size(dev{cond},2)
                    [a b] =  min(dev{cond}(:,pair));
                    [aa bb] = min(mean(devControl{cond}(:,pair,:),3));
                    imp = (a-mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(:,pair,:),3)));
                    zerolag = (min(dev{cond}(1:6,pair)) - mean(dev{cond}(:,pair))) ./ (aa - mean(mean(devControl{cond}(1,pair,:),3)));
                    if zerolag < 1 
                        zerolag = 1;
                    end
                    if imp > 4 & b > 7 & b < 150 &  zerolag < 1.1 & mean(dev{cond}(:,pair))>40
                       assemblyStrength = [assemblyStrength; imp];
                       optimalWindow = [optimalWindow;b];
%                        stren = [stren; imp imp2 b (min_mse_phase_all(pairs(pair,1),cond)-min_mse_rate(pairs(pair,1),cond)) ...
%                         ./ (min_mse_phase_all(pairs(pair,1),cond)+min_mse_rate(pairs(pair,1),cond)) ...
%                         (min_mse_phase_all(pairs(pair,1),cond))];
                        p = [p; pairs(pair,:)]; % cell pair passed threshold
                    elseif imp <= 4 & b > 7 & b < 150 &  zerolag < 1.1 & mean(dev{cond}(:,pair))>40
                        pp = [pp; pairs(pair,:)]; % cel pair did not pass threshold
                    end
                    pairCount = 1 + pairCount;
               end  
               
                f=[]; % get LS cell list
                for t=1:length(spikes.times)
                    if strcmp(spikes.region{t},'ls')
                    f = [f;t];
                    end
                end
                if ~isempty(p) % at least one assembly was found for this LS neuron
                    phaseRateBias = [phaseRateBias; (min_mse_phase_all(p(:,1),cond)-min_mse_rate(p(:,1),cond)) ];
                    ff = f(~ismember(f,p(:,1)));
                else
                if ~isempty(pp)
                    ff = pp(:,1);
                else
                    ff=[];
                end
                end
                    noAssembly_phaseRateBias = [noAssembly_phaseRateBias; (min_mse_phase_all(ff,cond)-min_mse_rate(ff,cond)) ];
               end
            end
        end
    end
    cd(homeDirectory)
end


%% plot summary data        
        
subplot(6,5,26)
scatter(optimalWindow,assemblyStrength,'.k')
axis([ 0 150 4 35])
hold on

subplot(6,5,27)
histogram(optimalWindow,1:2:200)
axis([1 1000 0 100])
set(gca,'xscale','log')

clear a n c ap np cp aa nn aap nnp
for iter = 1:1000
    r = randperm(length(phaseRateBias));
    aa(iter,:) = phaseRateBias(r(1:round(prctile(1:length(phaseRateBias),50))));
    r = randperm(length(noAssembly_phaseRateBias));
    nn(iter,:) = noAssembly_phaseRateBias(r(1:round(prctile(1:length(noAssembly_phaseRateBias),50))));
    for j = 1:100
    aap(iter,j) = prctile(aa(iter,:),j);
    nnp(iter,j) = prctile(nn(iter,:),j);
    end
end
subplot(6,5,29)
boundedline(1:100,mean(aap),3*std(aap),'m');
boundedline(1:100,mean(nnp),3*std(nnp),'k');
line([0 100],[0 0],'color','k')
axis([0 100 -6000 3000])



