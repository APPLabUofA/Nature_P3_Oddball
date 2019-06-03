%% Plot ERPs of your epochs
% An ERP is a plot of the raw, usually filtered, data, at one or multiple electrodes. It doesn't use time-frequency data.
% We make ERPs if we have segmented datasets that we want to compare across conditions.

% In this case, you are using all your electrodes, not a subset.
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
% % Pictures
% % erpdata(1,1) = Baseline
% % erpdata(2,1) = Nature
% % erpdata(3,1) = Urban

%%
%%%%%BELOW IS FOR PICTURE ONSET ALL PARTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = ['b','g','r'];
%%%%% Pick your condition %%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;

%%%%%ERPs for nature/urban pictures%%%%%
figure;hold on;
for i_part = 1:nparts
    subplot(ceil(sqrt(nparts)),ceil(sqrt(nparts)),i_part);boundedline(EEG.times,erpdata_parts(cond1).cond(:,i_part),std(erpdata_parts(cond1).cond(:,i_part))./sqrt(length(exp.participants)),col(cond1),...
        EEG.times,erpdata_parts(cond2).cond(:,i_part),std(erpdata_parts(cond2).cond(:,i_part))./sqrt(length(exp.participants)),col(cond2),...
        EEG.times,erpdata_parts(cond3).cond(:,i_part),std(erpdata_parts(cond3).cond(:,i_part))./sqrt(length(exp.participants)),col(cond3));
    
    xlim([-200 1000])
    ylim([-30 30])
    set(gca,'Ydir','reverse');
    line([0 0],[-10 10],'color','k');
    line([-200 1000],[0 0],'color','k');
    title(['Participant ' num2str(i_part)]);
end
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%BELOW IS FOR PICTURE ONSET GRAND AVERAGE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = ['b','g','r'];
%%%%% Pick your condition %%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;

%%%%%ERPs for nature/urban pictures%%%%%
figure;boundedline(EEG.times,erpdata(cond1).cond,std(erpdata_parts(cond1).cond,[],2)./sqrt(length(exp.participants)),col(cond1),...
    EEG.times,erpdata(cond2).cond,std(erpdata_parts(cond2).cond,[],2)./sqrt(length(exp.participants)),col(cond2),...
    EEG.times,erpdata(cond3).cond,std(erpdata_parts(cond3).cond,[],2)./sqrt(length(exp.participants)),col(cond3));

xlim([-200 1000])
ylim([-15 15])
set(gca,'Ydir','reverse');
line([0 0],[-10 10],'color','k');
line([-200 1000],[0 0],'color','k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%Power Topoplots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;

%%%%%Pick your time window 300-600 for P3, 150-250 for MMN%%%%%
% time1 = 150;
% time2 = 250;

time1 = 300;
time2 = 600;

%%%set limits for topoplots%%%
lims = [-2 2];
diff_lims = [-2 2];

%%%%%Get topoplots for differnce between targets and standards%%%%%
time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

erp_diff_out_bas = (all_chan_erpdata(cond1).cond);
erp_diff_out_nat = (all_chan_erpdata(cond2).cond);
erp_diff_out_urb = (all_chan_erpdata(cond3).cond);

figure('Color',[1 1 1]);
set(gca,'Color',[1 1 1]);
temp1 = mean(erp_diff_out_bas(:,time_window),2);
temp2 = mean(erp_diff_out_nat(:,time_window),2);
temp3 = mean(erp_diff_out_urb(:,time_window),2);
temp1(33:34) = NaN;
temp2(33:34) = NaN;
temp3(33:34) = NaN;

subplot(1,3,1);
topoplot((temp1),exp.electrode_locs, 'whitebk','on','plotrad',.6,'maplimits',lims)
title('Baseline');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');

subplot(1,3,2);
topoplot((temp2),exp.electrode_locs, 'whitebk','on','plotrad',.6,'maplimits',lims)
title('Nature');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');

subplot(1,3,3);
topoplot((temp3),exp.electrode_locs, 'whitebk','on','plotrad',.6,'maplimits',lims)
title('Urban');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');

%%%%%Difference Topoplots For Nature and Urbansud %%%%%
figure('Color',[1 1 1]);
set(gca,'Color',[1 1 1]);

topoplot((temp2-temp3),exp.electrode_locs, 'whitebk','on','plotrad',.6,'maplimits',diff_lims)
title('Nature - Urban');
t = colorbar('peer',gca);
set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Bar Plots For Time Windows%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;

%%%%%Pick your time window 300-600 for P3, 150-250 for MMN%%%%%
time1 = 300;
time2 = 600;

% % time1 = 150;
% % time2 = 250;

time_window = find(EEG.times>time1,1)-1:find(EEG.times>time2,1)-2;

baseline_pics = erpdata(cond1).cond(time_window);
nature_pics = erpdata(cond2).cond(time_window);
urban_pics = erpdata(cond3).cond(time_window);

%%%se across time
se_baseline = std(mean(erpdata_parts(cond1).cond(time_window,:),2),[],1)./sqrt(length(exp.participants));
se_nature = std(mean(erpdata_parts(cond2).cond(time_window,:),2),[],1)./sqrt(length(exp.participants));
se_urban = std(mean(erpdata_parts(cond3).cond(time_window,:),2),[],1)./sqrt(length(exp.participants));

%%%se across parts
% se_baseline = std(mean(erpdata_parts(cond1).cond(time_window,:),1),[],2)./sqrt(length(exp.participants));
% se_nature = std(mean(erpdata_parts(cond2).cond(time_window,:),1),[],2)./sqrt(length(exp.participants));
% se_urban = std(mean(erpdata_parts(cond3).cond(time_window,:),1),[],2)./sqrt(length(exp.participants));

figure; bar([mean(baseline_pics);mean(nature_pics);mean(urban_pics)]);
xlabel(['Baseline','Nature','Urban']);ylabel(['Mean Power from ' num2str(time1) ' to ' num2str(time2) ' at electrode ' electrodes(i_chan)]);

hold on;
errorbar([mean(baseline_pics);mean(nature_pics);mean(urban_pics)],...
    [se_baseline,se_nature,se_urban],'.');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%