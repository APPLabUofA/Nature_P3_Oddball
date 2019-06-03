electrodes = {EEG.chanlocs(:).labels};
% Type "electrodes" into the command line. This will show you which number to use for i_chan

% This code will take a grand average of the subjects, making one figure per set.
% This is the normal way to present ERP results.
i_chan = [3];%%%3 = Pz; 6 = Fz; 29 = F3; 30 = F4; 5 = FCz; 9 = PO3; 10 = PO4; 2 = Oz; 4 = Cz; 7 = O1; 8 = O2
for i_set = 1:nsets
    data_out = [];
    exp.setname{i_set}
    % The code below uses a nested loop to determine which segmented dataset corresponds to the right argument in data_out
    % e.g. if you have 5 sets, 20 participants, and 4 events, for i_set ==
    % 2 and i_part == 1 and i_event == 1, the code uses the data from set (2-1)*4*20 + (1-1)*20 + 1 == set 81
    for eegset = 1:nevents
        exp.event_names{1,eegset}
        for i_part = 1:nparts
            data_out(:,:,eegset,i_part,i_set) = nanmean(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).data,3);
            ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).filename
        end
        
        % this is the averaged ERP data. We take the mean across the selected channels (dim1) and the participants (dim4).
        erpdata(i_set,eegset).cond = squeeze(mean(mean(data_out(i_chan,:,eegset,:,i_set),1),4));
        erpdata_parts(i_set,eegset).cond = squeeze(mean(data_out(i_chan,:,eegset,:,i_set),1));
        all_chan_erpdata(i_set,eegset).cond = squeeze(mean(data_out(:,:,eegset,:,i_set),4));
        all_chan_erpdata_parts(i_set,eegset).cond = squeeze(data_out(:,:,eegset,:,i_set));
        
    end
end
%%%How data is organised
% % Low Tones
% % erpdata(1,1) = Baseline
% % erpdata(2,1) = Nature
% % erpdata(3,1) = Urban
% %
% % High Tones
% % erpdata(1,2) = Baseline
% % erpdata(2,2) = Nature
% % erpdata(3,2) = Urban


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%T-Tests tone onset during specified window%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;
%%%%%Pick your tone 1 = low 2 = high%%%%%
tone1 = 1;
tone2 = 2;
%%%%%Pick your time window 350-550 for P3, 150-250 for MMN%%%%%
time_window = 2;

if time_window == 1
    
    time1 = 150;
    time2 = 250;
    
elseif time_window == 2
    time1 = 350;
    time2 = 550;
    
end

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

%%%Calculate difference wave for two conditions%%%
%%%Need to calculate the average voltage across your time window, for each
%%%participant. Will be comparing list of participant means.
erp_diff_1 = (erpdata_parts(cond1,tone2).cond-erpdata_parts(cond1,tone1).cond);
erp_diff_2 = (erpdata_parts(cond2,tone2).cond-erpdata_parts(cond2,tone1).cond);

%%%left tailed t test%%%
[h_l1 p_l1 ci_l1 stat_l1] = ttest(mean(erp_diff_1(time_window,:),1),0,.05,'left',2)
[h_l2 p_l2 ci_l2 stat_l2] = ttest(mean(erp_diff_2(time_window,:),1),0,.05,'left',2)

%%%right tailed t test%%%
[h_r1 p_r1 ci_r1 stat_r1] = ttest(mean(erp_diff_1(time_window,:),1),0,.05,'right',2)
[h_r2 p_r2 ci_r2 stat_r2] = ttest(mean(erp_diff_2(time_window,:),1),0,.05,'right',2)

%%%two tailed comparing nature and urban%%%
[h_b p_b ci_b stat_b] = ttest(mean(erp_diff_1(time_window,:),1),mean(erp_diff_2(time_window,:),1),.05,'both',2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%T-Tests picture onset at specified window%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 1;
cond2 = 2;
%%%%%Pick your time window 300-600, 150-250%%%%%
time_window = 1;

if time_window == 1
    
    time1 = 150;
    time2 = 250;
    
elseif time_window == 2
    
    time1 = 300;
    time2 = 600;
    
end

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

%%%left tailed t test%%%
[h_l1 p_l1 ci_l1 stat_l1] = ttest(mean(erpdata_parts(cond1).cond(time_window,:),1),0,.05,'left',2)
[h_l2 p_l2 ci_l2 stat_l2] = ttest(mean(erpdata_parts(cond2).cond(time_window,:),1),0,.05,'left',2)

%%%right tailed t test%%%
[h_r1 p_r1 ci_r1 stat_r1] = ttest(mean(erpdata_parts(cond1).cond(time_window,:),1),0,.05,'right',2)
[h_r2 p_r2 ci_r2 stat_r2] = ttest(mean(erpdata_parts(cond2).cond(time_window,:),1),0,.05,'right',2)

%%%two tailed t test comparing conditions%%%
[h_b p_b ci_b stat_b] = ttest(mean(erpdata_parts(cond1).cond(time_window,:),1),mean(erpdata_parts(cond2).cond(time_window,:),1),.05,'both',2)

%%
%%%grab trial numbers%%%

baseline_low = [];
baseline_high = [];
nature_low = [];
nature_high = [];
urban_low = [];
urban_high = [];

count1 = 1;
count2 = 1;
count3 = 1;
count4 = 1;
count5 = 1;
count6 = 1;

for i_event = 1:length(ALLEEG)
    
    if strcmp(ALLEEG(i_event).filename(end-21:end),'Baseline_Low_Tones.set') == 1
        baseline_low(count1) = length(ALLEEG(i_event).event);
        count1 = count1+1;
    elseif strcmp(ALLEEG(i_event).filename(end-22:end),'Baseline_High_Tones.set') == 1
        baseline_high(count2) = length(ALLEEG(i_event).event);
        count2 = count2+1;
    elseif strcmp(ALLEEG(i_event).filename(end-19:end),'Nature_Low_Tones.set') == 1
        nature_low(count3) = length(ALLEEG(i_event).event);
        count3 = count3+1;
    elseif strcmp(ALLEEG(i_event).filename(end-20:end),'Nature_High_Tones.set') == 1
        nature_high(count4) = length(ALLEEG(i_event).event);
        count4 = count4+1;
    elseif strcmp(ALLEEG(i_event).filename(end-18:end),'Urban_Low_Tones.set') == 1
        urban_low(count5) = length(ALLEEG(i_event).event);
        count5 = count5+1;
    elseif strcmp(ALLEEG(i_event).filename(end-19:end),'Urban_High_Tones.set') == 1
        urban_high(count6) = length(ALLEEG(i_event).event);
        count6 = count6+1;
    end
end


mean(baseline_low)
mean(baseline_high)
mean(nature_low)
mean(nature_high)
mean(urban_low)
mean(urban_high)

min(baseline_low)
min(baseline_high)
min(nature_low)
min(nature_high)
min(urban_low)
min(urban_high)

max(baseline_low)
max(baseline_high)
max(nature_low)
max(nature_high)
max(urban_low)
max(urban_high)

%%
%%%%%STATS FOR REVIEWS%%%%%
%%%times for P3 and MMN/P2 windows%%%
i_cond = 1; %1 for tones; 2 for pics

if i_cond == 1
    time1_p3 = 350;
    time2_p3 = 550;
    
    time1_mmn = 150;
    time2_mmn = 250;
    
elseif i_cond == 2
    
    time1_p3 = 300;
    time2_p3 = 600;
    
    time1_mmn = 100;
    time2_mmn = 250;
    
end

time_window_p3 = find(EEG.times>time1_p3,1)-1:find(EEG.times>time2_p3,1)-2;
time_window_mmn = find(EEG.times>time1_mmn,1)-1:find(EEG.times>time2_mmn,1)-2;
for i_chan = 1:length(electrodes)
    
    if i_cond == 1
        baseline_lowtones = squeeze(all_chan_erpdata_parts(1,1).cond(i_chan,:,:));
        nature_lowtones = squeeze(all_chan_erpdata_parts(2,1).cond(i_chan,:,:));
        urban_lowtones = squeeze(all_chan_erpdata_parts(3,1).cond(i_chan,:,:));
        
        baseline_hightones = squeeze(all_chan_erpdata_parts(1,2).cond(i_chan,:,:));
        nature_hightones = squeeze(all_chan_erpdata_parts(2,2).cond(i_chan,:,:));
        urban_hightones = squeeze(all_chan_erpdata_parts(3,2).cond(i_chan,:,:));
        
        baseline_diff = baseline_hightones - baseline_lowtones;
        nature_diff = nature_hightones - nature_lowtones;
        urban_diff = urban_hightones - urban_lowtones;
        
    elseif i_cond == 2
        
        baseline_diff = squeeze(all_chan_erpdata_parts(1,1).cond(i_chan,:,:));
        nature_diff = squeeze(all_chan_erpdata_parts(2,1).cond(i_chan,:,:));
        urban_diff = squeeze(all_chan_erpdata_parts(3,1).cond(i_chan,:,:));
        
    end
    
    %%%Effect Size Comparing Conditions%%%
    %%%P3%%%
    stats_bn_p3(i_chan) = mes([mean(baseline_diff(time_window_p3,:),1)]',[mean(nature_diff(time_window_p3,:),1)]','hedgesg','isDep',1);
    stats_bu_p3(i_chan) = mes([mean(baseline_diff(time_window_p3,:),1)]',[mean(urban_diff(time_window_p3,:),1)]','hedgesg','isDep',1);
    stats_nu_p3(i_chan) = mes([mean(nature_diff(time_window_p3,:),1)]',[mean(urban_diff(time_window_p3,:),1)]','hedgesg','isDep',1);
    
    %%%MMN%%%
    stats_bn_mmn(i_chan) = mes([mean(baseline_diff(time_window_mmn,:),1)]',[mean(nature_diff(time_window_mmn,:),1)]','hedgesg','isDep',1);
    stats_bu_mmn(i_chan) = mes([mean(baseline_diff(time_window_mmn,:),1)]',[mean(urban_diff(time_window_mmn,:),1)]','hedgesg','isDep',1);
    stats_nu_mmn(i_chan) = mes([mean(nature_diff(time_window_mmn,:),1)]',[mean(urban_diff(time_window_mmn,:),1)]','hedgesg','isDep',1);
    
    %%%Effect Size for Individual Conditions%%%
    %%%P3%%%
    stats_b_p3(i_chan) = mes([mean(baseline_diff(time_window_p3,:),1)]',[0],'hedgesg');
    stats_n_p3(i_chan) = mes([mean(nature_diff(time_window_p3,:),1)]',[0],'hedgesg');
    stats_u_p3(i_chan) = mes([mean(urban_diff(time_window_p3,:),1)]',[0],'hedgesg');
    
    %%%MMN%%%
    stats_b_mmn(i_chan) = mes([mean(baseline_diff(time_window_mmn,:),1)]',[0],'hedgesg');
    stats_n_mmn(i_chan) = mes([mean(nature_diff(time_window_mmn,:),1)]',[0],'hedgesg');
    stats_u_mmn(i_chan) = mes([mean(urban_diff(time_window_mmn,:),1)]',[0],'hedgesg');
    
    %%%Power Comparing Conditions%%%
    %%%P3%%%
    pwrout_bn_p3(i_chan,1) = sampsizepwr('t2',[mean(mean(baseline_diff(time_window_p3,:),1)), std(mean(baseline_diff(time_window_p3,:),1))],...
        mean(mean(nature_diff(time_window_p3,:),1)),[],length(mean(baseline_diff(time_window_p3,:),1)));
    
    pwrout_bu_p3(i_chan,1) = sampsizepwr('t2',[mean(mean(baseline_diff(time_window_p3,:),1)), std(mean(baseline_diff(time_window_p3,:),1))],...
        mean(mean(urban_diff(time_window_p3,:),1)),[],length(mean(baseline_diff(time_window_p3,:),1)));
    
    pwrout_nu_p3(i_chan,1) = sampsizepwr('t2',[mean(mean(nature_diff(time_window_p3,:),1)), std(mean(nature_diff(time_window_p3,:),1))],...
        mean(mean(urban_diff(time_window_p3,:),1)),[],length(mean(baseline_diff(time_window_p3,:),1)));
    
    %%%MMN%%%
    pwrout_bn_mmn(i_chan,1) = sampsizepwr('t2',[mean(mean(baseline_diff(time_window_mmn,:),1)), std(mean(baseline_diff(time_window_mmn,:),1))],...
        mean(mean(nature_diff(time_window_mmn,:),1)),[],length(mean(baseline_diff(time_window_mmn,:),1)));
    
    pwrout_bu_mmn(i_chan,1) = sampsizepwr('t2',[mean(mean(baseline_diff(time_window_mmn,:),1)), std(mean(baseline_diff(time_window_mmn,:),1))],...
        mean(mean(urban_diff(time_window_mmn,:),1)),[],length(mean(baseline_diff(time_window_mmn,:),1)));
    
    pwrout_nu_mmn(i_chan,1) = sampsizepwr('t2',[mean(mean(nature_diff(time_window_mmn,:),1)), std(mean(nature_diff(time_window_mmn,:),1))],...
        mean(mean(urban_diff(time_window_mmn,:),1)),[],length(mean(baseline_diff(time_window_mmn,:),1)));
    
    %%%Power For Individual Conditions%%%
    %%%P3%%%
    pwrout_b_p3(i_chan,1) = sampsizepwr('t',[mean(mean(baseline_diff(time_window_p3,:),1)), std(mean(baseline_diff(time_window_p3,:),1))],...
        0,[],length(mean(baseline_diff(time_window_p3,:),1)));
    
    pwrout_n_p3(i_chan,1) = sampsizepwr('t',[mean(mean(nature_diff(time_window_p3,:),1)), std(mean(nature_diff(time_window_p3,:),1))],...
        0,[],length(mean(nature_diff(time_window_p3,:),1)));
    
    pwrout_u_p3(i_chan,1) = sampsizepwr('t',[mean(mean(urban_diff(time_window_p3,:),1)), std(mean(urban_diff(time_window_p3,:),1))],...
        0,[],length(mean(urban_diff(time_window_p3,:),1)));
    
    %%%MMN%%%
    pwrout_b_mmn(i_chan,1) = sampsizepwr('t',[mean(mean(baseline_diff(time_window_mmn,:),1)), std(mean(baseline_diff(time_window_mmn,:),1))],...
        0,[],length(mean(baseline_diff(time_window_mmn,:),1)));
    
    pwrout_n_mmn(i_chan,1) = sampsizepwr('t',[mean(mean(nature_diff(time_window_mmn,:),1)), std(mean(nature_diff(time_window_mmn,:),1))],...
        0,[],length(mean(nature_diff(time_window_mmn,:),1)));
    
    pwrout_u_mmn(i_chan,1) = sampsizepwr('t',[mean(mean(urban_diff(time_window_mmn,:),1)), std(mean(urban_diff(time_window_mmn,:),1))],...
        0,[],length(mean(urban_diff(time_window_mmn,:),1)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Perform T-tests for entire ERP duration--TONE ONSET%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;
%%%%%Pick your tone 1 = low 2 = high%%%%%
tone1 = 1;
tone2 = 2;

col = ['b','b';'g','g';'r','r'];

ttest_times = zeros(length(EEG.times),4);
erp_diff_1 = (erpdata_parts(cond2,tone2).cond-erpdata_parts(cond2,tone1).cond);
erp_diff_2 = (erpdata_parts(cond3,tone2).cond-erpdata_parts(cond3,tone1).cond);

for i_time = 1:length(EEG.times)
    
    [h p ci stat] = ttest(mean(erp_diff_1(i_time,:),1),mean(erp_diff_2(i_time,:)),.05,'both',2);
    
    ttest_times(i_time,1)= h;
    ttest_times(i_time,2)= p;
    
end

ttest_sig = ttest_times(:,2);
ttest_non = ttest_times(:,2);

ttest_sig(ttest_sig > 0.05) = NaN;

ttest_times(ttest_times == 0 ) = NaN;

figure;hold on;

boundedline(EEG.times(800:2000),(erpdata(cond1,tone2).cond(800:2000)-erpdata(cond1,tone1).cond(800:2000)),std(erpdata_parts(cond1,tone2).cond(800:2000,:)-erpdata_parts(cond1,tone1).cond(800:2000,:),[],2)./sqrt(length(exp.participants)),col(cond1,tone1),...
    EEG.times(800:2000),(erpdata(cond2,tone2).cond(800:2000)-erpdata(cond2,tone1).cond(800:2000)),std(erpdata_parts(cond2,tone2).cond(800:2000,:)-erpdata_parts(cond2,tone1).cond(800:2000,:),[],2)./sqrt(length(exp.participants)),col(cond2,tone2),...
    EEG.times(800:2000),(erpdata(cond3,tone2).cond(800:2000)-erpdata(cond3,tone1).cond(800:2000)),std(erpdata_parts(cond3,tone2).cond(800:2000,:)-erpdata_parts(cond3,tone1).cond(800:2000,:),[],2)./sqrt(length(exp.participants)),col(cond3,tone2));

plot(EEG.times(800:2000),(ttest_times(800:2000,1)*14),'m','LineWidth',6);

xlim([-200 1000])
ylim([-15  15])
set(gca,'Ydir','reverse');
line([0 0],[-10 10],'color','k');
line([-200 1000],[0 0],'color','k');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TTEST FOR TONE ONSET AT EACH ELECTRODE, AT SPECIFIED TIME WINDOWS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;
%%%%%Pick your tone 1 = low 2 = high%%%%%
tone1 = 1;
tone2 = 2;
%%%%%Pick your time window 350-550 for P3, 150-250 for MMN%%%%%

time_window = 1;

if time_window == 1
    
    time1 = 150;
    time2 = 250;
    
elseif time_window == 2
    
    time1 = 350;
    time2 = 550;
    
end

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

ttest_electrodes = zeros(length(electrodes),2);
ttest_electrodes_ci = zeros(length(electrodes),2);
ttest_electrodes_stats = [];
ttest_electrodes_stats.electrode = zeros(length(electrodes),1);

for i_elec = 1:length(electrodes)
    
    %%%Calculate difference wave for two conditions%%%
    %%%Need to calculate the average voltage across your time window, for each
    %%%participant. Will be comparing list of participant means.
    
    erp_diff_1 = (squeeze(all_chan_erpdata_parts(cond1,tone2).cond(i_elec,:,:))-squeeze(all_chan_erpdata_parts(cond1,tone1).cond(i_elec,:,:)));
    erp_diff_2 = (squeeze(all_chan_erpdata_parts(cond2,tone2).cond(i_elec,:,:))-squeeze(all_chan_erpdata_parts(cond2,tone1).cond(i_elec,:,:)));
    
    [h p ci stat] = ttest(mean(erp_diff_1(time_window,:),1),mean(erp_diff_2(time_window,:),1),.05,'both',2)
    
    ttest_electrodes(i_elec,1) = h;
    ttest_electrodes(i_elec,2) = p;
    ttest_electrodes_ci(i_elec,1) = ci(1);
    ttest_electrodes_ci(i_elec,2) = ci(2);
    ttest_electrodes_stats(i_elec,1).electrode = stat;
    
end
sigchans = find(ttest_electrodes(:,1) == 1);

%%%%%Get topoplots for differnce between targets and standards%%%%%

erp_diff_out_nat = (all_chan_erpdata(cond1,tone2).cond-all_chan_erpdata(cond1,tone1).cond);
erp_diff_out_urb = (all_chan_erpdata(cond2,tone2).cond-all_chan_erpdata(cond2,tone1).cond);

figure('Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
temp1 = mean(erp_diff_out_nat(:,time_window),2);
temp2 = mean(erp_diff_out_urb(:,time_window),2);
temp1(33:34) = NaN;
temp2(33:34) = NaN;

%%%%%Difference Topoplots%%%%%
topoplot((temp1-temp2),exp.electrode_locs, 'whitebk','on','plotrad',.6,'maplimits',[-2 2],'emarker2',{sigchans,'*','w'});
title('Nature - Urban');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TTEST FOR PICTURE ONSET AT EACH ELECTRODE, AT SPECIFIED TIME WINDOWS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_chan = [2];

electrodes = {EEG.chanlocs(:).labels};
for i_chan = 1:length(electrodes)
    for i_set = 1:nsets
        data_out = [];
        exp.setname{i_set}
        % The code below uses a nested loop to determine which segmented dataset corresponds to the right argument in data_out
        % e.g. if you have 5 sets, 20 participants, and 4 events, for i_set ==
        % 2 and i_part == 1 and i_event == 1, the code uses the data from set (2-1)*4*20 + (1-1)*20 + 1 == set 81
        for eegset = 1:nevents
            exp.event_names{1,eegset}
            for i_part = 1:nparts
                % % % % %             data_out(i_chan,time,set,part,event)
                % % % % %             ALLEEG(1).data(i_chan,time,trial)
                %%%%%Look at the entire 3 blocks for each condition%%%%%
                data_out(:,:,eegset,i_part,i_set) = nanmean(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).data,3);
                
                ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).filename
            end
            
            % this is the averaged ERP data. We take the mean across the selected channels (dim1) and the participants (dim4).
            erpdata(i_chan,i_set,eegset).cond = squeeze(mean(mean(data_out(i_chan,:,eegset,:,i_set),1),4));
            erpdata_parts(i_chan,i_set,eegset).cond = squeeze(mean(data_out(i_chan,:,eegset,:,i_set),1));
            all_chan_erpdata(i_set,eegset).cond = squeeze(mean(data_out(:,:,eegset,:,i_set),4));
            %             all_chan_erpdata(i_chan,i_set,eegset).cond = squeeze(mean(mean(data_out(i_chan,:,eegset,:,i_set),1),4));
            %     erpdata = squeeze(mean(mean(data_out(i_chan,:,:,:,i_set),1),4));
            
        end
    end
end


%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;

%%%%%Pick your time window 100-250, 300-600%%%%%
time1 = 100;
time2 = 250;

time1 = 300;
time2 = 600;

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

ttest_electrodes = zeros(length(electrodes),2);
ttest_electrodes_adjusted = zeros(length(electrodes),2);
ttest_electrodes_p_corrected_adjusted = zeros(length(electrodes),2);
ttest_electrodes_ci = zeros(length(electrodes),2);
ttest_electrodes_stats.electrode = zeros(length(electrodes),1);

for i_elec = 1:length(electrodes)
    
    % %     [h p ci stat] = ttest(mean(erpdata_parts(i_elec,cond1).cond(time_window,:),1),mean(erpdata_parts(i_elec,cond2).cond(time_window,:),1),.025,'both',2);
    [h p ci stat] = ttest(mean(erpdata_parts(i_elec,cond1).cond(time_window,:),1),mean(erpdata_parts(i_elec,cond2).cond(time_window,:),1),.05,'both',2);
    
    ttest_electrodes(i_elec,1) = h;
    ttest_electrodes(i_elec,2) = p;
    ttest_electrodes_ci(i_elec,1) = ci(1);
    ttest_electrodes_ci(i_elec,2) = ci(2);
    ttest_electrodes_stats(i_elec,1).electrode = stat;
    
end

sigchans = find(ttest_electrodes(:,1) == 1);

%%%%%Get topoplots for differnce between targets and standards%%%%%

erp_diff_out_nat = (all_chan_erpdata(cond1).cond);
erp_diff_out_urb = (all_chan_erpdata(cond2).cond);

figure('Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
temp1 = mean(erp_diff_out_nat(:,time_window),2);
temp2 = mean(erp_diff_out_urb(:,time_window),2);
temp1(33:34) = NaN;
temp2(33:34) = NaN;

%%%%%Difference Topoplots%%%%%
topoplot((temp1-temp2),exp.electrode_locs, 'whitebk','on','plotrad',.6,'maplimits',[-2 2],'emarker2',{sigchans,'*','w'});
title('Nature - Urban');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TTEST FOR TONE ONSET AT EACH ELECTRODE, ACROSS ENTIRE ERP%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% i_chan = [3];
% 
% electrodes = {EEG.chanlocs(:).labels};
% for i_chan = 1:length(electrodes)
%     for i_set = 1:nsets
%         data_out = [];
%         exp.setname{i_set}
%         % The code below uses a nested loop to determine which segmented dataset corresponds to the right argument in data_out
%         % e.g. if you have 5 sets, 20 participants, and 4 events, for i_set ==
%         % 2 and i_part == 1 and i_event == 1, the code uses the data from set (2-1)*4*20 + (1-1)*20 + 1 == set 81
%         for eegset = 1:nevents
%             exp.event_names{1,eegset}
%             for i_part = 1:nparts
%                 
%                 data_out(:,:,eegset,i_part,i_set) = nanmean(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).data,3);
%                 
%                 ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).filename
%             end
%             
%             % this is the averaged ERP data. We take the mean across the selected channels (dim1) and the participants (dim4).
%             erpdata(i_chan,i_set,eegset).cond = squeeze(mean(mean(data_out(i_chan,:,eegset,:,i_set),1),4));
%             erpdata_parts(i_chan,i_set,eegset).cond = squeeze(mean(data_out(i_chan,:,eegset,:,i_set),1));
%             all_chan_erpdata(i_set,eegset).cond = squeeze(mean(data_out(:,:,eegset,:,i_set),4));
%             
%         end
%     end
% end


%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;
%%%%%Pick your tone 1 = low 2 = high%%%%%
tone1 = 1;
tone2 = 2;

%%%need to store tests for each electrode, and at each time point%%%
ttest_electrodes_h = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_h_corrected = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_p = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_p_corrected = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_ci_lower = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_ci_upper = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_stats = [];
ttest_electrodes_stats.electrode = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_t_stats = zeros(length(electrodes),length(EEG.times));

for i_elec = 1:length(electrodes)
    for i_time = 1:length(EEG.times)
        
        %%%For this we are going to compare the difference between nature
        %%%and urban difference waves, and then compare those values to
        %%%zero
        %%%Calculate difference wave for nature and urban, and then calculate the difference of that%%%
        %%%Need to calculate the average voltage across your time window, for each
        %%%participant. Will be comparing list of participant means.
        %%%nature difference wave%%%
%         erp_diff_1 = (erpdata_parts(i_elec,cond1,tone2).cond-erpdata_parts(i_elec,cond1,tone1).cond);
        erp_diff_1 = (squeeze(all_chan_erpdata_parts(cond1,tone2).cond(i_elec,:,:))-squeeze(all_chan_erpdata_parts(cond1,tone1).cond(i_elec,:,:)));
        %%%urban difference wave%%%
%         erp_diff_2 = (erpdata_parts(i_elec,cond2,tone2).cond-erpdata_parts(i_elec,cond2,tone1).cond);
        erp_diff_2 = (squeeze(all_chan_erpdata_parts(cond2,tone2).cond(i_elec,:,:))-squeeze(all_chan_erpdata_parts(cond2,tone1).cond(i_elec,:,:)));
        %%%Difference between nature and urban differene waves%%%
%         erp_diff_3 = (erpdata_parts(i_elec,cond1,tone2).cond-erpdata_parts(i_elec,cond1,tone1).cond)-(erpdata_parts(i_elec,cond2,tone2).cond-erpdata_parts(i_elec,cond2,tone1).cond);
        erp_diff_3 = (squeeze(all_chan_erpdata_parts(cond1,tone2).cond(i_elec,:,:))-squeeze(all_chan_erpdata_parts(cond1,tone1).cond(i_elec,:,:)))-...
            (squeeze(all_chan_erpdata_parts(cond2,tone2).cond(i_elec,:,:))-squeeze(all_chan_erpdata_parts(cond2,tone1).cond(i_elec,:,:)));
        
        [h p ci stat] = ttest(erp_diff_3(i_time,:),0,.05,'both',2);
        % %     [h p ci stat] = ttest(mean(erp_diff_1(time_window,:),1),mean(erp_diff_2(time_window,:),1),.025,'right',2);
        
        ttest_electrodes_h(i_elec,i_time) = h;
        ttest_electrodes_p(i_elec,i_time) = p;
        ttest_electrodes_ci_lower(i_elec,i_time) = ci(1);
        ttest_electrodes_ci_upper(i_elec,i_time) = ci(2);
        ttest_electrodes_stats(i_elec,i_time).electrode = stat;
        ttest_electrodes_t_stats(i_elec,i_time) = stat.tstat;
    end
    [ttest_electrodes_p_corrected(i_elec,:), h_corrected] = bonf_holm(ttest_electrodes_p(i_elec,:), 0.05);
    ttest_electrodes_h_corrected(i_elec,:) = h_corrected;
end

colormap redblue
sig_tone = ttest_electrodes_t_stats;

cutoff = tinv(.025,length(exp.participants)-1);%for a two tailed test with alpha=0.05
clus_lims_tone = find(diff(sum(abs(sig_tone)>abs(cutoff))>3));

figure;imagesc(EEG.times(800:2000),2:32,ttest_electrodes_h(2:32,(800:2000)));colorbar;set(gca,'YTick',[2:32]); colormap redblue;
figure;imagesc(EEG.times(800:2000),2:32,(1-ttest_electrodes_p(2:32,(800:2000))));colorbar;set(gca,'YTick',[2:32]); colormap redblue;

%%%this is a graph showing all the ttest results over the entire ERP
%%%following tone onset%%%
figure;hold on;
imagesc(EEG.times(800:2000),2:32,sig_tone(2:32,(800:2000)));colormap(b2r(cutoff*5,abs(cutoff)*5));colorbar;set(gca,'YTick',[2:32],'YTickLabel',electrodes(2:32));
ylabel('Electrode');xlabel('Time (ms)'); title('T-test of Difference Waves (Nature - Urban)');

%%%same graph as above but only showing the significant times%%%
figure;hold on;
sig_tone(ttest_electrodes_t_stats>cutoff&ttest_electrodes_t_stats<abs(cutoff))= 0;
imagesc(EEG.times(800:2000),2:32,sig_tone(2:32,(800:2000)));colormap(b2r(cutoff,abs(cutoff)));colorbar;set(gca,'YTick',[2:32],'YTickLabel',electrodes(2:32));
ylabel('Electrode');xlabel('Time (ms)'); title('T-test of Difference Waves (Nature - Urban)'); set(gca, 'YDir','reverse'); hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TTEST FOR PICTURE ONSET AT EACH ELECTRODE, ACROSS ENTIRE ERP%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_chan = [3];

electrodes = {EEG.chanlocs(:).labels};
for i_chan = 1:length(electrodes)
    for i_set = 1:nsets
        data_out = [];
        exp.setname{i_set}
        % The code below uses a nested loop to determine which segmented dataset corresponds to the right argument in data_out
        % e.g. if you have 5 sets, 20 participants, and 4 events, for i_set ==
        % 2 and i_part == 1 and i_event == 1, the code uses the data from set (2-1)*4*20 + (1-1)*20 + 1 == set 81
        for eegset = 1:nevents
            exp.event_names{1,eegset}
            for i_part = 1:nparts
                % % % % %             data_out(i_chan,time,set,part,event)
                % % % % %             ALLEEG(1).data(i_chan,time,trial)
                %%%%%Look at the entire 3 blocks for each condition%%%%%
                data_out(:,:,eegset,i_part,i_set) = nanmean(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).data,3);
                
                ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).filename
            end
            
            % this is the averaged ERP data. We take the mean across the selected channels (dim1) and the participants (dim4).
            erpdata(i_chan,i_set,eegset).cond = squeeze(mean(mean(data_out(i_chan,:,eegset,:,i_set),1),4));
            erpdata_parts(i_chan,i_set,eegset).cond = squeeze(mean(data_out(i_chan,:,eegset,:,i_set),1));
            all_chan_erpdata(i_set,eegset).cond = squeeze(mean(data_out(:,:,eegset,:,i_set),4));
            %             all_chan_erpdata(i_chan,i_set,eegset).cond = squeeze(mean(mean(data_out(i_chan,:,eegset,:,i_set),1),4));
            %     erpdata = squeeze(mean(mean(data_out(i_chan,:,:,:,i_set),1),4));
            
        end
    end
end


%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;

%%%need to store tests for each electrode, and at each time point%%%
ttest_electrodes_h = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_p = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_ci_lower = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_ci_upper = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_stats.electrode = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_t_stats = zeros(length(electrodes),length(EEG.times));

for i_elec = 1:length(electrodes)
    for i_time = 1:length(EEG.times)
        
        erp_diff = erpdata_parts(i_elec,cond1).cond - erpdata_parts(i_elec,cond2).cond;
        
        % %     [h p ci stat] = ttest(mean(erpdata_parts(i_elec,cond1).cond(time_window,:),1),mean(erpdata_parts(i_elec,cond2).cond(time_window,:),1),.025,'both',2);
        [h p ci stat] = ttest(erp_diff(i_time,:),0,.05,'both',2);
        
        ttest_electrodes_h(i_elec,i_time) = h;
        ttest_electrodes_p(i_elec,i_time) = p;
        ttest_electrodes_ci_lower(i_elec,i_time) = ci(1);
        ttest_electrodes_ci_upper(i_elec,i_time) = ci(2);
        ttest_electrodes_stats(i_elec,i_time).electrode = stat;
        ttest_electrodes_t_stats(i_elec,i_time) = stat.tstat;
        
    end
end
sig_pic = ttest_electrodes_t_stats;

cutoff = tinv(.025,length(exp.participants)-1);%for a two tailed test with alpha=0.05
clus_lims_pic = find(diff(sum(abs(sig_pic)>abs(cutoff))>3));

%%%this is a graph showing all the ttest results over the entire ERP
%%%following tone onset%%%
figure;imagesc(EEG.times(800:2000),2:32,sig_pic(2:32,(800:2000)));colormap(b2r(cutoff*5,abs(cutoff)*5));colorbar;set(gca,'YTick',[2:32],'YTickLabel',electrodes(2:32));
ylabel('Electrode');xlabel('Time (ms)'); title('T-test of Difference Waves (Nature - Urban)');line([0,0],[0,32]);
for ii = 1:length(clus_lims_pic); line([clus_lims_pic(ii)+.5 clus_lims_pic(ii)+.5],[0 33],'color','k'); end

%%%same graph as above but only showing the significant times%%%
sig_pic(ttest_electrodes_t_stats>cutoff&ttest_electrodes_t_stats<abs(cutoff))=0;
figure;imagesc(EEG.times(800:2000),2:32,sig_pic(2:32,(800:2000)));colormap(b2r(cutoff,abs(cutoff)));colorbar;set(gca,'YTick',[2:32],'YTickLabel',electrodes(2:32));
ylabel('Electrode');xlabel('Time (ms)'); title('T-test of Difference Waves (Nature - Urban)');line([0,0],[0,32]);
for ii = 1:length(clus_lims_pic); line([clus_lims_pic(ii)+.5 clus_lims_pic(ii)+.5],[0 33],'color','k'); end
