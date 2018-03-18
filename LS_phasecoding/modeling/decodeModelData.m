clear all
% close all

key = {'dynamic_wo_theta_seq_ls_output_','dynamic_wo_theta_seq_trial_';...
       'dynamic_ls_output_','dynamic_trial_';...
       'random_uniform_ls_output_','random_uniform_trial_';...
       'static_ls_output_','static_trial_';};
   
for fff = 1:4
    figure(fff)
    clf
    ls_train_files = dir(['data/' key{fff,1} '*']);
    % theta_files = dir('data/*theta*');
    % [b a] = butter(3,[5/500 12/500],'bandpass');

    for trial = 1:length(ls_train_files)
        dat = load(['data/' key{fff,1}  num2str(trial-1) '.mat']);
        ls_trains{trial} = dat.block.segments{end}.spiketrains{1}.times;
        dat = load(['data/' key{fff,2}  num2str(trial-1) '_theta.mat']);
        theta(trial,:) = makeLength(dat.theta,6500);
        angles(trial,:) = angle(hilbert((zscore(theta(trial,:)))));
    end

    smoothing = 20*21;
    positionDecodingMaxCorr_binned_box_mean.results{1} = table;


    for smoothing = 21
    cond = 1;
    nBins = 200;
    cell=1;

    phase_trains_smooth=[];
    cos_phase_trains_smooth=[];
    sin_phase_trains_smooth=[];
    rates_trains_smooth = [];
    pos=[];

    for t = 1:length(ls_train_files)
        rates = zeros(200,1);
        phases = zeros(200,1);
        rates(round(ls_trains{t}./32.5))=1;
        f = round(ls_trains{t});
        ff = round(ls_trains{t}./32.5);
        phases(ff)=angles(t,f);
        % cut first XXX ms off of simulation....
        phases = circ_smoothTS(phases,smoothing,'method','mean','exclude',0);
        rates = smooth(rates,smoothing);    
        rates = (rates(25:175))';
        phases = (phases(25:175))';

        phase_trains_smooth=[phase_trains_smooth;...
            makeLength(phases,200)'];

        rates_trains_smooth = [rates_trains_smooth; ...
                               makeLength(rates,200)'];    
        pos = [pos;[1:nBins]'];                   
    end 

                r = rates_trains_smooth;
                p = phase_trains_smooth;
                p_cos =cos(phase_trains_smooth);
                p_sin = sin(phase_trains_smooth);

                count = 1;
                for iter = 1:10
                rr = randperm(length(r));
                pct = round(prctile(1:length(r),50));

                r_train = r(rr(1:pct));
                r_test = r(rr(pct+1:end));
                p_train = p(rr(1:pct));
                p_test = p(rr(pct+1:end));
                p_cos_train = p_cos(rr(1:pct));
                p_cos_test = p_cos(rr(pct+1:end));
                p_sin_train = p_sin(rr(1:pct));
                p_sin_test = p_sin(rr(pct+1:end));


                pos_train = pos(rr(1:pct));
                pos_test = pos(rr(pct+1:end));
                % rate 
                cl = max_correlation_coefficient_CL;
                cl = train(cl,[r_train'],pos_train');
                yfit = test(cl,[r_test']);
                struct.mse_rate = nanmean((pos_test'-yfit).^2);
    %             struct.mse_rate_pval = stats.p';
                % phase all
                cl = max_correlation_coefficient_CL;
                cl = train(cl,[p_cos_train'; p_sin_train';p_train'],pos_train');
                yfit = test(cl,[p_cos_test'; p_sin_test';p_test']);
                struct.mse_phase_all = nanmean((pos_test'-yfit).^2);
    %             struct.mse_phase_all_pval = stats.p';
                %phase
                cl = max_correlation_coefficient_CL;
                cl = train(cl,[p_train'],pos_train');
                yfit = test(cl,[p_test']);
                struct.mse_phase = nanmean((pos_test'-yfit).^2);
    %             struct.mse_phase_pval = stats.p';
                %phase cos
                cl = max_correlation_coefficient_CL;
                cl = train(cl,[p_cos_train'],pos_train');
                yfit = test(cl,[p_cos_test']);
                struct.mse_phase_cos = nanmean((pos_test'-yfit).^2);
    %             struct.mse_phase_cos_pval = stats.p';
                %phase sin
    %             cl = max_correlation_coefficient_CL;
    %             cl = train(cl,[p_sin_train'],pos_train');
    %             yfit = test(cl,[p_sin_test']);
    %             struct.mse_phase_sin = nanmedian((pos_test'-yfit).^2);
    %             struct.mse_phase_sin_pval = stats.p';
                % all predictors
                cl = max_correlation_coefficient_CL;
                cl = train(cl,[p_cos_train'; p_sin_train';p_train';r_train'],pos_train');
                yfit = test(cl,[p_cos_test'; p_sin_test';p_test';r_test']);
                struct.mse_both  = nanmean((pos_test'-yfit).^2);
    %             struct.mse_both_pval = stats.p';

                cl = max_correlation_coefficient_CL;
                cl = train(cl,[r_train'],pos_train');
                yfit = test(cl,[r_test(randperm(length(r_test)))']);
                struct.mse_chance_rate  = nanmean((pos_test'-yfit).^2);

                cl = max_correlation_coefficient_CL;
                cl = train(cl,[p_cos_train'; p_sin_train';p_train'],pos_train');
                yfit = test(cl,[p_cos_test(randperm(length(p_test)))'; p_sin_test(randperm(length(p_test)))';p_test(randperm(length(p_test)))']);
                struct.mse_chance_phase  = nanmean((pos_test'-yfit).^2);

                % store peripherals
                struct.iter = count;
    %             struct.dfe = stats.dfe;
                struct.tau = smoothing;
                struct.condition = cond;
                positionDecodingMaxCorr_binned_box_mean.results{cell} = [positionDecodingMaxCorr_binned_box_mean.results{cell}; struct2table(struct)];

                count = 1+count;
                end
    %             if cell == 205 && cond == 2
    % % %                 
    % %                 figure(cond)
                    t_rate = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_rate',...
                    'GroupingVariables',{'tau','condition'});
                    t_phase_all = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase_all',...
                    'GroupingVariables',{'tau','condition'});
                    t_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase',...
                    'GroupingVariables',{'tau','condition'});
                    t_phase_cos = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_phase_cos',...
                    'GroupingVariables',{'tau','condition'});
                    t_both = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_both',...
                    'GroupingVariables',{'tau','condition'});
                    t_chance_rate = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_rate',...
                    'GroupingVariables',{'tau','condition'});   
                    t_chance_phase = varfun(@mean,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_phase',...
                    'GroupingVariables',{'tau','condition'});  
                    t_std_rate = varfun(@std,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_rate',...
                    'GroupingVariables',{'tau','condition'});  
                    t_std_phase = varfun(@std,positionDecodingMaxCorr_binned_box_mean.results{cell},'InputVariables','mse_chance_phase',...
                    'GroupingVariables',{'tau','condition'});  
                    tab = join(join(join(join(join(join(t_rate,t_phase_all),t_both),t_chance_rate),t_phase),t_phase_cos),t_chance_phase);
                    rows = find(tab.condition==cond);
                    subplot(3,2,1);
                    plot(tab.tau(rows),tab.mean_mse_phase(rows),'.g')
                    hold on
                    plot(tab.tau(rows),tab.mean_mse_phase_cos(rows),'g')

    %                 imagesc(phase_trains_smooth);
    %                 caxis([-pi pi])
                    subplot(3,2,2);
                    boundedline(tab.tau(rows),tab.mean_mse_chance_rate(rows),3*t_std_rate.std_mse_chance_rate(rows),'k')
                    hold on
    %                 nans(smoothing) = sum(isnan(yfit));
                    plot(tab.tau(rows),tab.mean_mse_rate(rows),'r')
                    plot(tab.tau(rows),tab.mean_mse_both(rows),'k')
                    plot(tab.tau(rows),tab.mean_mse_phase(rows),'g')
                    subplot(3,2,3)
                    boundedline(tab.tau(rows),tab.mean_mse_chance_phase(rows),3*t_std_phase.std_mse_chance_phase(rows),'k')
                    hold on

    %                 plot(positionDecodingMaxCorr_binned_box_mean.results{cell}.tau(rows),...
    %                     mean(positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_rate(rows,:)'),'r')
    %                 hold on
    %                 plot(positionDecodingMaxCorr_binned_box_mean.results{cell}.tau(rows),...
    %                     mean(positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_phase_all(rows,:)'),'g')
                    plot(positionDecodingMaxCorr_binned_box_mean.results{cell}.tau(rows),...
                        mean(positionDecodingMaxCorr_binned_box_mean.results{cell}.mse_both(rows,:)'),'k')

    %                 subplot(3,2,3)
    %                 imagesc(rates_trains_smooth);
    %                 subplot(3,2,4)

                    title([cell cond smoothing]')

                    pause(.1)



                    end
    subplot(3,2,4)

    p = positionDecodingMaxCorr_binned_box_mean;
    rows = find(p.results{1}.tau==21);
    rc = p.results{1}.mse_chance_rate(rows);
    r = p.results{1}.mse_rate(rows);
    pp = p.results{1}.mse_phase(rows);
    [a rrc] = kstest2(r,rc);
    [a ppc] = kstest2(pp,rc);
    [a ppr] = kstest2(pp,r);


    scatter(ones(10,1)*2,p.results{1}.mse_rate(rows),'.r')
    hold on
    scatter(ones(10,1)*3,p.results{1}.mse_phase(rows),'.b')
    scatter(ones(10,1),p.results{1}.mse_chance_rate(rows),'.k')
    axis([0 4 0 9000])
    errorbar(1,mean(p.results{1}.mse_chance_rate(rows)),std(p.results{1}.mse_chance_rate(rows)),'k')
    errorbar(2,mean(p.results{1}.mse_rate(rows)),std(p.results{1}.mse_rate(rows)),'r')
    errorbar(3,mean(p.results{1}.mse_phase(rows)),std(p.results{1}.mse_phase(rows)),'b')
    sigstar({[1,2], [1,3], [2,3]},[rrc ppc ppr])

    subplot(3,2,5)

    smoothing = 21;
    rates_trains_smooth = [];
    for t = 1:length(ls_train_files)
        rates = zeros(200,1);

        rates(round(ls_trains{t}./32.5))=1;

        rates = smooth(rates,smoothing);    
        rates = (rates(25:175))';

        rates_trains_smooth = [rates_trains_smooth, ...
                               makeLength(rates,200)'];    

    end 
%     imagesc(rates_trains_smooth')
plot(mean(rates_trains_smooth'))
    title(key{fff,1})

end
