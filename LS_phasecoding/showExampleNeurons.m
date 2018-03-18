% script for looking through other cool neurons

% recording, cell #, behavioral condition, region

examples = {'20170505_396um_0um_merge',                          80,5,'hpc';
            '20170501_216um_0um_170501_112210',                  44,7,'ls';
            '20170509_468um_36um_170509_103451',                 2, 2,'ls';
            '20170509_468um_36um_170509_103451',                 24,1,'ls';
            '20170504_396um_0um_merge',                          19,1,'ls';
            '20170525_900um_936um_170525_150124' ,               8, 2,'ls'
              };
%           
% examples = {'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,1,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,2,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,3,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,4,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,5,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,6,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,7,'ls';
% 'DT2_rPPC_rCCG_3540um_1288um_20160227_160227_121226',44,8,'ls';
% };
%           %             
% %           412 cell 205 cond 2
% % 525  cell 7 cond 2
%           % linear track examples
% otherExamples = {'20170502_288um_0um_merge',                 4, 2,'ls';
%                      '20170507_468um_36um_merge',                3, 1,'ls';
%                      '20170507_468um_36um_merge',                18,1,'ls';
%                      '20170507_468um_36um_merge',                3, 2,'ls';
%                      '20170507_468um_36um_merge',                8, 2,'ls';
%                      '20170508_468um_36um_170508_102251',        28,2,'ls';
%                      '20170522_900um_936um_170522_151132',       4, 1,'ls';
%                      '20170522_900um_936um_170522_151132',       12, 2,'ls'};
%            % circular track examples
%            
% examples = {'20161021_1840um_2412um_merge'   ,         [  9]  ,  [1]  ,  'ls'
% '20161021_1840um_2412um_merge'    ,        [  9]  ,  [2]  ,  'ls'
% '20161021_1840um_2412um_merge'    ,        [ 41]  ,  [1]  ,  'ls'
% '20161021_1840um_2412um_merge'    ,        [ 44]  ,  [2]  ,  'ls'
% '20161021_1840um_2412um_merge'     ,       [ 52]  ,  [1]  ,  'ls'
% '20161026_1944um_2232um_merge'     ,       [  4]  ,  [2]  ,  'ls'
% '20161026_1944um_2232um_merge'      ,      [  6]  ,  [1]  ,  'ls'
% '20161026_1944um_2232um_merge'      ,      [  6]  ,  [2]  ,  'ls'
% '20161026_1944um_2232um_merge'      ,      [ 12]  ,  [1]  ,  'ls'
% '20161026_1944um_2232um_merge'      ,      [ 21]  ,  [2]  ,  'ls'
% '20161026_1944um_2232um_merge'      ,      [ 25]  ,  [1]  ,  'ls'
% '20161026_1944um_2232um_merge'      ,      [ 25]  ,  [2]  ,  'ls'
% '20170412_1440um_1170um_170412_103153' ,   [205]  ,  [2]  ,  'ls'
% '20170502_288um_0um_merge'           ,     [ 37]  ,  [1]  ,  'ls'
% '20170504_396um_0um_merge'           ,     [  5]  ,  [2]  ,  'ls'
% '20170504_396um_0um_merge'           ,     [ 22]  ,  [1]  ,  'ls'
% '20170504_396um_0um_merge'           ,     [ 29]  ,  [2]  ,  'ls'
% '20170521_900um_936um_170521_125602' ,     [ 16] ,   [1] ,   'ls'
% '20170522_900um_936um_170522_151132' ,     [  3]  ,  [1]  ,  'ls'
% '20170522_900um_936um_170522_151132' ,     [ 11] ,  [2],   'ls'
% '20170522_900um_936um_170522_151132' ,    [ 12]  ,  [1]  ,  'ls'
% '20170522_900um_936um_170522_151132' ,     [ 14] ,   [1] ,   'ls'
% '20170522_900um_936um_170522_151132' ,     [ 14],[2]    ,'ls'
% '20170522_900um_936um_170522_151132' ,     [ 21] ,   [2] ,   'ls'
% '20170522_900um_936um_170522_151132'      [ 23]  ,  [1]  ,  'ls'
% '20170525_900um_936um_170525_150124' ,     [  7] ,   [2] ,   'ls'
% '20170525_900um_936um_170525_150124' ,     [  8] ,   [1] ,   'ls'
% '20170525_900um_936um_170525_150124' ,     [  8] ,   [2] ,   'ls'
% '20170525_900um_936um_170525_150124' ,     [  9] ,   [1] ,   'ls'
% '20170525_900um_936um_170525_150124' ,     [ 16] ,   [2] ,   'ls'};


% for ex = 1:size(examples,1)
%     cd(examples{ex,1})
d = dir('*201*');
for i = 45:length(d)
    cd(d(i).name)
    % load necessary data
    sessionInfo = bz_getSessionInfo;
    load([sessionInfo.FileName '.behavior.mat'])
    load([sessionInfo.FileName '.firingMaps.cellinfo.mat'])
    load([sessionInfo.FileName '.phaseMaps.cellinfo.mat'])
    load([sessionInfo.FileName '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
%     load([examples{ex,1} '.behavior.mat'])
%     load([examples{ex,1} '.firingMaps.cellinfo.mat'])
%     load([examples{ex,1} '.phaseMaps.cellinfo.mat'])
%     load([examples{ex,1} '.positionDecodingMaxCorr_binned_box_mean.cellinfo.mat'])
    
    spikes = bz_GetSpikes; % requires buzcode
    for cell=1:length(spikes.times)
        for cond = 1:length(unique(behavior.events.trialConditions))
            if sum(behavior.events.trialConditions==cond)>10
%     cond = examples{ex,3};
%     cell = examples{ex,2};
    
    % start plotting
%     figure(ex)
    clf
    subplot(3,3,1)
    bz_plotTrials(behavior,'condition',cond)
%     scatter(phaseMaps.phaseMaps{cond}{cell}(:,3),phaseMaps.phaseMaps{cond}{cell}(:,4),'.m')
    axis square
    ylabel(['position (' behavior.units ')'])
    ylabel(['position (' behavior.units ')'])
    
    subplot(3,3,2)
    imagesc(squeeze(firingMaps.rateMaps{cond}(cell,:,:)))
    ylabel('trial #')
    xlabel('linearized position')
    
    subplot(3,3,3)
    if ~isempty(phaseMaps.phaseMaps{cond}{cell}) % prevents [] empty error
    scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),(phaseMaps.phaseMaps{cond}{cell}(:,end)),'.k');
    hold on
    scatter(phaseMaps.phaseMaps{cond}{cell}(:,1),(phaseMaps.phaseMaps{cond}{cell}(:,end))+2*pi,'.k');
    nBins = length(behavior.events.map{1}.x);
    axis([0 nBins -pi pi*3])
    xlabel('linearized position')
    ylabel('theta phase')
    end
    
    subplot(3,3,4)
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
    [a rc] = kstest2(r,c);
    [a rp] = kstest2(r,p);
    [a pc] = kstest2(p,c);
    sigstar({[2,3],[1,2], [1,3]},[rp rc pc])
%     title(['recording: ' examples{ex,1}])
    
    subplot(3,3,5)
    plot(spikes.rawWaveform{cell}*.195)
    ylabel('millivolts')
    xlabel('sample #')
    
    subplot(3,3,6)
    histogram(diff(spikes.times{cell}),0:.001:.03)
    ylabel('count')
    xlabel('time (seconds)')
    
    
    
    
%     saveFigure(['/home/david/Dropbox/Documents/pubs/ls_precession/' examples{ex,1} '_' num2str(cell)])
    

if pc < .001 & rc > .01 & length(phaseMaps.phaseMaps{cond}{cell}) > 40 & strcmp(phaseMaps.region{cell},'ls')
%     pause
    savefig(['/home/david/Dropbox/temp/' sessionInfo.FileName '_' num2str(cell) '_' num2str(cond) '.fig'])
    clf
end
% savefig(['/home/david/Dropbox/Documents/pubs/ls_precession/matlab_figs/1/' examples{ex,1} '_' num2str(examples{ex,2}) '.fig']) 
% cd D:\Dropbox\datasets\lsDataset
    
%     pause
clear p
            end
        end
    end
    cd /home/david/datasets/lsDataset/
end