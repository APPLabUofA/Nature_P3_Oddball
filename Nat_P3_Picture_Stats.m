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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Perform T-tests for entire ERP duration--PICTURE ONSET%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;

col = ['b','b';'g','g';'r','r'];

ttest_times = zeros(length(EEG.times),4);

for i_time = 1:length(EEG.times)
    [h p ci stat] = ttest(mean(erpdata_parts(cond2).cond(i_time,:),1),mean(erpdata_parts(cond3).cond(i_time,:),1),.05,'both',2);
    
    ttest_times(i_time,1)= h;
    ttest_times(i_time,2)= p;
    
end

ttest_sig = ttest_times(:,2);
ttest_non = ttest_times(:,2);

ttest_sig(ttest_sig > 0.05) = NaN;

ttest_times(ttest_times == 0 ) = NaN;

figure;hold on;
boundedline(EEG.times(800:2000),erpdata(cond1).cond(800:2000),std(erpdata_parts(cond1).cond(800:2000,:),[],2)./sqrt(length(exp.participants)),col(cond1),...
    EEG.times(800:2000),erpdata(cond2).cond(800:2000),std(erpdata_parts(cond2).cond(800:2000,:),[],2)./sqrt(length(exp.participants)),col(cond2),...
    EEG.times(800:2000),erpdata(cond3).cond(800:2000),std(erpdata_parts(cond3).cond(800:2000,:),[],2)./sqrt(length(exp.participants)),col(cond3));

plot(EEG.times(800:2000),(ttest_times(800:2000,1)*14),'m','LineWidth',6);

xlim([-200 1000])
ylim([-15  15])
set(gca,'Ydir','reverse');
%%%%% epoched to last entrainer %%%%%
line([0 0],[-10 10],'color','k');
line([-200 1000],[0 0],'color','k');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%TTEST FOR PICTURE ONSET AT EACH ELECTRODE, AT SPECIFIED TIME WINDOWS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;

%%%%%Pick your time window 100-250, 300-600%%%%%
time1 = 100;
time2 = 250;

time1 = 300;
time2 = 600;

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

clear ttest_electrodes_h ttest_electrodes_p ttest_electrodes_ci_lower ttest_electrodes_ci_upper ttest_electrodes_stats ttest_electrodes_t_stats

ttest_electrodes = zeros(length(electrodes),2);
ttest_electrodes_adjusted = zeros(length(electrodes),2);
ttest_electrodes_p_corrected_adjusted = zeros(length(electrodes),2);
ttest_electrodes_ci = zeros(length(electrodes),2);
ttest_electrodes_stats.electrode = zeros(length(electrodes),1);

for i_elec = 1:length(electrodes)
    
    [h p ci stat] = ttest(mean(squeeze(all_chan_erpdata_parts(cond1).cond(i_elec,time_window,:)),1),mean(squeeze(all_chan_erpdata_parts(cond2).cond(i_elec,time_window,:)),1),.05,'both',2);
    
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
%%%TTEST FOR PICTURE ONSET AT EACH ELECTRODE, ACROSS ENTIRE ERP%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 2;
cond2 = 3;

%%%need to store tests for each electrode, and at each time point%%%
clear ttest_electrodes_h ttest_electrodes_p ttest_electrodes_ci_lower ttest_electrodes_ci_upper ttest_electrodes_stats ttest_electrodes_t_stats

ttest_electrodes_h = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_p = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_ci_lower = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_ci_upper = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_stats.electrode = zeros(length(electrodes),length(EEG.times));
ttest_electrodes_t_stats = zeros(length(electrodes),length(EEG.times));

for i_elec = 1:length(electrodes)
    for i_time = 1:length(EEG.times)
        
        erp_diff = squeeze(all_chan_erpdata_parts(cond1).cond(i_elec,:,:)) - squeeze(all_chan_erpdata_parts(cond2).cond(i_elec,:,:));
        
        [h p ci stat] = ttest(erp_diff(i_time,:),0,.025,'both',2);
        
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
ylabel('Electrode');xlabel('Time (ms)'); title('T-test of Difference Waves (Nature - Urban)');

%%%same graph as above but only showing the significant times%%%
sig_pic(ttest_electrodes_t_stats>cutoff&ttest_electrodes_t_stats<abs(cutoff))=0;
figure;imagesc(EEG.times(800:2000),2:32,sig_pic(2:32,(800:2000)));colormap(b2r(cutoff,abs(cutoff)));colorbar;set(gca,'YTick',[2:32],'YTickLabel',electrodes(2:32));
ylabel('Electrode');xlabel('Time (ms)'); title('T-test of Difference Waves (Nature - Urban)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

