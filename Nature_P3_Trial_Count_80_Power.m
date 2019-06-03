%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Calculate Trial Count for 80% Power%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pick_trials_stand = [ 5 10:10:150 ]; %150 for tones; 200 for visual
pick_trials_targ = pick_trials_stand/5;
pick_trials_stand = pick_trials_stand-pick_trials_targ;

stim_type = 2; %1 for tones, 2 for pics
time_window = 1; %ONLY FOR VISUAL; 1 = 100-250ms; 2 = 300-600ms
test_type = 1; %1 = left-tailed t-test; 2 = right-tailed t-test
electrode = 1;% 3 = pz; 6 = fz; 2 = oz
perms = 10000;

if stim_type == 1
    if electrode == 6
        window = [350 550];
    elseif electrode == 3
        window = [150 250];
    end
elseif stim_type == 2
    if electrode == 1
        if time_window == 1
            window = [100 250];
        elseif time_window == 2
            window = [300 600];
        end
    elseif electrode == 6
        if time_window == 1
            window = [100 250];
        elseif time_window == 2
            window = [300 600];
        end
    end
end

all_stand_trials = [];
all_targ_trials = [];
for i_part = 1:nparts
    for i_set = 1:nsets
        exp.setname{i_set}
        for eegset = 1:nevents
            exp.event_names{1,eegset}
            ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).filename
            
            if i_part == 1 && i_set == 1 && eegset == 1
                time_window = find(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).times>=window(1),1):find(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).times>=window(2),1)-1; %data points of time window
                erp_out = zeros(length(time_window),2,nparts,nsets,perms,length(pick_trials_stand));
                
            end
            
            if eegset == 1
                n_stand_trials = ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).trials;
                all_stand_trials(i_part,i_set) = n_stand_trials;
                for i_perm = 1:perms
                    for i_pick = 1:length(pick_trials_stand)
                        
                        stand_trials = randperm(n_stand_trials,pick_trials_stand(i_pick));            %without replacement (permutation)
                        %                 stand_trials = randi(n_stand_trials,pick_trials_stand(i_pick));  %with relacement (bootstrap)
                        erp_out(:,1,i_part,i_set,i_perm,i_pick) = mean(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).data(electrode,time_window,stand_trials),3);
                    end
                end
                 
            elseif eegset == 2
                n_targ_trials = ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).trials;
                all_targ_trials(i_part,i_set) = n_targ_trials;
                for i_perm = 1:perms
                    for i_pick = 1:length(pick_trials_stand)
                        
                        targ_trials = randperm(n_targ_trials,pick_trials_targ(i_pick));
                        %                   targ_trials = randi(n_targ_trials,pick_trials_targ(i_pick),1)';
                        
                        erp_out(:,2,i_part,i_set,i_perm,i_pick) = mean(ALLEEG((i_set-1)*nevents*nparts + (eegset-1)*(nparts) + i_part).data(electrode,time_window,targ_trials),3);
                    end
                    
                end
            end
            
        end
        
    end
end
%%%erp_out(times ,events ,parts, conditions, perms, picks)
%%%mean_erp(events ,parts, conditions, perms, picks)
%%%mean_erp_difft(parts, conditions, perms, picks)

H_base =[];
P_base = [];
H_nat =[];
P_nat = [];
H_urb =[];
P_urb = [];
mean_erp = squeeze(mean(erp_out,1));
mean_erp_diff = squeeze(mean_erp(2,:,:,:,:) - mean_erp(1,:,:,:,:));
for i_pick = 1:length(pick_trials_stand)
    for i_perm = 1:perms
        
        if test_type == 1
            [H_base(1,i_pick,i_perm),P_base(1,i_pick,i_perm)] =...
                ttest(mean_erp_diff(:,1,i_perm,i_pick),0,.05,'left',1);%%%for MMN
            [H_nat(1,i_pick,i_perm),P_nat(1,i_pick,i_perm)] =...
                ttest(mean_erp_diff(:,2,i_perm,i_pick),0,.05,'left',1);%%%for MMN
            [H_urb(1,i_pick,i_perm),P_urb(1,i_pick,i_perm)] =...
                ttest(mean_erp_diff(:,3,i_perm,i_pick),0,.05,'left',1);%%%for MMN
            
        elseif test_type == 2
            
            [H_base(1,i_pick,i_perm),P_base(1,i_pick,i_perm)] = ttest(mean_erp_diff(:,1,i_perm,i_pick),0,.05,'right',1);%%%for P3
            [H_nat(1,i_pick,i_perm),P_nat(1,i_pick,i_perm)] = ttest(mean_erp_diff(:,2,i_perm,i_pick),0,.05,'right',1);%%%for P3
            [H_urb(1,i_pick,i_perm),P_urb(1,i_pick,i_perm)] = ttest(mean_erp_diff(:,3,i_perm,i_pick),0,.05,'right',1);%%%for P3
            %
        end
    end
end
prop_sig_base = sum(H_base,3)/perms;
prop_sig_nat = sum(H_nat,3)/perms;
prop_sig_urb = sum(H_urb,3)/perms;


figure; hold on;
plot(pick_trials_targ(end:-1:1),prop_sig_base(:,end:-1:1), 'color', 'b');
plot(pick_trials_targ(end:-1:1),prop_sig_nat(:,end:-1:1), 'color', 'g');
plot(pick_trials_targ(end:-1:1),prop_sig_urb(:,end:-1:1), 'color', 'r');
legend('Baseline','Nature','Urban','location','SouthEast');  line([0 135],[.8 .8],'color','k'); axis tight;
plot(pick_trials_targ(end:-1:1), sqrt(pick_trials_targ(end:-1:1))/max(sqrt(pick_trials_targ(end:-1:1))),'k');
xlim([0,40]); hold off;

figure; hold on;
plot(pick_trials_stand(end:-1:1),prop_sig_base(:,end:-1:1),'-bo');
plot(pick_trials_stand(end:-1:1),prop_sig_nat(:,end:-1:1),'-go');
plot(pick_trials_stand(end:-1:1),prop_sig_urb(:,end:-1:1),'-ro');
legend('Baseline','Nature','Urban','location','SouthEast');  line([0 400],[.8 .8],'color','k'); axis tight;
plot(pick_trials_stand(end):-1:1, sqrt(pick_trials_stand(end):-1:1)/max(sqrt(pick_trials_stand(end):-1:1)),'k'); hold off;
xlim([0,150]); hold off;
