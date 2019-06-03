%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Perform T-tests for entire ERP duration--PICTURE ONSET%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Pick your condition 1 = baseline 2 = nature 3 = urban%%%%%
cond1 = 1;
cond2 = 2;
cond3 = 3;

col = ['b','b';'g','g';'r','r'];
% % col = ['g','g';'r','r'];

ttest_times = zeros(length(EEG.times),4);

for i_time = 1:length(EEG.times)
    [h p ci stat] = ttest(mean(erpdata_parts(cond2).cond(i_time,:),1),mean(erpdata_parts(cond3).cond(i_time,:),1),.05,'both',2);
    % % [h p ci stat] = ttest(mean(erpdata_parts(cond1).cond(i_time,:),1),mean(erpdata_parts(cond2).cond(i_time,:),1),.025,'both',2);
    
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

% plot(EEG.times,ttest_non,'k');
% plot(EEG.times,ttest_sig,'m');
plot(EEG.times(800:2000),(ttest_times(800:2000,1)*14),'m','LineWidth',6);

xlim([-200 1000])
ylim([-15  15])
set(gca,'Ydir','reverse');
%%%%% epoched to last entrainer %%%%%
line([0 0],[-10 10],'color','k');
line([-200 1000],[0 0],'color','k');
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%