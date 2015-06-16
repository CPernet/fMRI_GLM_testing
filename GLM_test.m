% Data analysis workflow to compare multiple GLM under two different fMRI
% conditions.
% Voxel 1 - differences in the amplitude of the neural response between 2 conditions.
% Voxel 2 - differences in the duration of the neural response between 2 conditions
% Each time check differences at the subject level in terms of data fir
% (R2) and the group level (average t value and power)
%
% Cyril Pernet May 2015

MC = 10000;
MeanR2 = NaN(3,9,MC); % average fit over 20 subjects, voxel 1/2 for 9 GLM repeated 1000 times (ie 2*20000 fits)
MeanDvalue  = NaN(3,13,MC); % group mean beta difference value cond1<cond2, voxel 1/2 from 9 GLM (+4 corrections) repeated 1000 times (ie 1000 gp effect)
Power_detection = NaN(3,13,MC); % binary result (significant/not signifiant) for the group analaysis from 9 GLM (+4 corrections) repeated 1000 times
Power_difference = NaN(3,13,MC); % binary result (significant/not signifiant) for the group analaysis from 9 GLM (+4 corrections) repeated 1000 times

for simulation =1:MC
    fprintf('MC %g\n',simulation)
    
    % 1 - generate RT
    RT(simulation) = OHBM2015_generateRT;
    
    % 2 - fit voxel 1
    [MeanR2(1,:,simulation), MeanDvalue(1,:,simulation), Power_detection(1,:,simulation),Power_difference(1,:,simulation)] = OHBM2015_intensity_voxel(RT(simulation));
    
    % 3 - fit voxel 2
    [MeanR2(2,:,simulation), MeanDvalue(2,:,simulation), Power_detection(2,:,simulation),Power_difference(2,:,simulation)]= OHBM2015_duration_voxel(RT(simulation));

    % 4 - fit voxel 3
    [MeanR2(3,:,simulation), MeanDvalue(3,:,simulation), Power_detection(3,:,simulation),Power_difference(3,:,simulation)]= OHBM2015_fixed_voxel(RT(simulation));

end
save all 

% mean R2 per model
figure; 
violin(squeeze(MeanR2(1,[1 3 5 7 2 4 6 8 9],:))'); title('Mean R2, amplitude changes trial to trial')
axis tight; grid on
figure; 
violin(squeeze(MeanR2(2,[1 3 5 7 2 4 6 8 9],:))'); title('Mean R2, duration changes trial to trial')
axis tight; grid on
figure; 
violin(squeeze(MeanR2(3,[1 3 5 7 2 4 6 8 9],:))'); title('Mean R2, no trial to trial change')
axis tight; grid on

% Power per model

tmp = nansum(Power_difference,3); 
tmp = tmp(:,[1 10 3 11 5 12 7 13 9]);

figure
[temp, CI] = binofit(tmp(3,:),10000);
bar(temp); title('Power when there is no trial to trial change'); 
hold on; errorbar([1:9],temp,temp-CI(:,1)',CI(:,2)'-temp,'rx','LineWidth',2);
grid on; axis([0.5 9.5 0 1.01])

figure
[temp, CI] = binofit(tmp(2,:),10000);
bar(temp); title('Power when there are duration changes'); 
hold on; errorbar([1:9],temp,temp-CI(:,1)',CI(:,2)'-temp,'rx','LineWidth',2);
grid on; axis([0.5 9.5 0 1.01])

figure
[temp, CI] = binofit(tmp(1,:),10000);
bar(temp); title('Power when there are amplitude changes'); 
hold on; errorbar([1:9],temp,temp-CI(:,1)',CI(:,2)'-temp,'rx','LineWidth',2);
grid on; axis([0.5 9.5 0 1.01])

% power with BF corrected only
figure
bar(mean(Power_difference([3 1 2],[1 10 3 11 5 12 7 13 9],:),3))
title('Power to detect a difference'); grid on

figure
bar(mean(Power_detection([3 1 2],[1 10 3 11 5 12 7 13 9],:),3))
title('Power to detect a activations'); grid on

figure
subplot(1,3,1); bar(mean(mean(Power_detection([3 1 2],[1 10 3 11 5 12 7 13 9],:),3)))
subplot(1,3,2); bar(mean(mean(Power_difference([3 1 2],[1 10 3 11 5 12 7 13 9],:),3)))
subplot(1,3,3); bar(mean((mean(Power_detection([3 1 2],[1 10 3 11 5 12 7 13 9],:),3)+mean(Power_difference([3 1 2],[1 10 3 11 5 12 7 13 9],:),3))./2))


MeanDvalue = squeeze(MeanDvalue(1,1,:));
[MeanDvalue,sorting_index]=sort(MeanDvalue);
bins = MeanDvalue(1:MC/10:MC);

figure; hold on
for model = 1:13
    tmp= squeeze(Power_detection(3,model,:)); tmp = tmp(sorting_index)';
    index = 1; index2 = 10;
    for b = 1:10;
        tmp2(b) = mean(tmp(index:index2));
        index = index+99; index2=index2+99;
    end
    subplot(4,4,model);plot(bins,tmp2)
end

% so what is the best model? given we don't know if we'll have more voxel
% that do show changes, or amplitude or shape changes - a safe option is to
% take the model which, on average, is the most powerfull for both
% detection above 0 and for testing differences.
