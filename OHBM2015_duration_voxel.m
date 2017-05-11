function [MeanR2, RT_diff, Power_detection, Power_difference]=OHBM2015_duration_voxel(RT)

% this simulation is designed to explore the effect
% of variable GLM on a voxel where neural amplitude varies with RT
% model 1: hrf
% model 2: BF
% model 3: hrf + PM 
% model 4: BF + PM 
% model 5: hrf mean RT
% model 6: BF mean RT
% model 7: hrf mean RT + PM
% model 8: BF mean RT + PM
% model 9: hrf with durations = RT
% --------------------------------
% Cyril Pernet 09-06-2014

% hemodynamic response function
% hrf model using SPM function
% ---------------------------
xBF.dt = 0.5;
xBF.name = 'hrf' (with time and dispersion derivative)';
xBF.length = 32; % over a 20 sec window
xBF.order = 1;
xBF.T = 30;
xBF2 = spm_get_bf(xBF);
xBF.UNITS = 'secs';


for subject =1:20
    RT_condition1 = RT(1,:,subject);
    RT_condition2 = RT(2,:,subject);
    RT_diff(subject) = mean(RT_condition2 - RT_condition1);
    
    %% BOLD simulation
    % now make 2 voxels with 2 conditions
    % voxel 1 duration = 0
    % voxel 2 duration = RT
    % note that for ease of interpretation
    % I work using the super-sampled design
    
    % make different BOLD data
    onsets = [1 31 81 93 161 201 218 291 321 361];
    duration = round(mean([RT_condition1 RT_condition2])*2); % * 2 because xBF.dt = 0.5
    
    Y1 = zeros(500,1); % condition 1, with 0ms duration
    Y2 = zeros(500,1); % condition 2, with 0ms duration
    X1 = zeros(500,1); % hrf (repeat twice)
    X2 = zeros(500,1); % hrf mean duration (repeat twice)
    X3 = zeros(500,1); % hrf duration = RT condition 1
    X4 = zeros(500,1); % hrf duration = RT condition 2
    for i=1:10
        X1(onsets(i)) = 1; % design matrix with duration 0
        X2(onsets(i):(onsets(i)+duration-1)) = 1; % design matrix with duration = mean RT
        stop = onsets(i) + round(RT_condition1(i)*2) -1;
        Y1(onsets(i):stop) = 1; % epoch = RT
        X3(onsets(i):stop) = 1; % design matrix for cond1 duration = RT
        stop = onsets(i) + round(RT_condition2(i)*2) -1;
        Y2(onsets(i):stop) = 1; % epoch = RT
        X4(onsets(i):stop) = 1; % design matrix for cond2 duration = RT
    end
    
    Y1 = conv(Y1+randn(500,1)./10,xBF.bf(:,1)); Y1 = Y1(1:400)+100;
    Y2 = conv(Y2+randn(500,1)./10,xBF.bf(:,1)); Y2 = Y2(1:400)+100;
    % mean(Y2-Y1)

    %% data modelling
    % For simplicity we put the 2 conditions one after the other Y1/Y2
    % we fit X to each one. It the conditions are orthogonal, then this
    % makes no difference from being inter-mixed ; of course estimates
    % differ if there is a correlation.
    
    % data
    Y = [[Y1 ; Y2]+randn(800,1)./2];
    SStotal = norm(Y-mean(Y)).^2;
    % [~,~,~,stat]=ttest(Y2-Y1);
    
    % model 1: hrf
    x1 = conv(X1,xBF.bf(:,1)); x1 = x1(1:400)+100;
    X = [[x1;zeros(400,1)] [zeros(400,1);x1] ones(800,1)];
    b = pinv(X)*Y; Beta1(subject,:) = b(1:2);
    Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R21(subject) = SSeffect / SStotal;
    
    % model 2: BF
    x11 = conv(X1,xBF.bf(:,1)); x11 = x11(1:400)+100;
    x12 = conv(X1,xBF.bf(:,2)); x12 = x12(1:400)+100;
    x13 = conv(X1,xBF.bf(:,3)); x13 = x13(1:400)+100;
    X = [[[x11 x12 x13];zeros(400,3)] [zeros(400,3);[x11 x12 x13]] ones(800,1)];
    X = [[[x11 x13];zeros(400,2)] [zeros(400,2);[x11 x13]] ones(800,1)];
    b = pinv(X)*Y; Beta2(subject,:) = b([1 4]);
    b1corrected = (sqrt(((b(1)^2)*sum(X(:,1).^2))+((b(2)^2)*sum(X(:,2).^2))+((b(3)^2)*sum(X(:,3).^2))))*sign(b(1));
    b4corrected = (sqrt(((b(4)^2)*sum(X(:,4).^2))+((b(5)^2)*sum(X(:,5).^2))+((b(6)^2)*sum(X(:,6).^2))))*sign(b(4));
    Beta2corrected(subject,:) = [b1corrected b4corrected]; Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R22(subject) = SSeffect / SStotal;
    
    % model 3: hrf + PM (a priori wining model)
    x1 = conv(X1,xBF.bf(:,1)); x1 = x1(1:400)+100;
    x2 = X1; x2(X1==1) = detrend(RT_condition1,'constant'); % parametric regressor cond1
    x2 = conv(x2,xBF.bf(:,1)); x2 = x2(1:400)+100;
    XPM1 = spm_orth([x1 x2]);
    x2 = X1; x2(X1==1) = detrend(RT_condition2,'constant'); % parametric regressor cond2
    x2 = conv(x2,xBF.bf(:,1)); x2 = x2(1:400)+100;
    XPM2 = spm_orth([x1 x2]);
    X = [[XPM1;zeros(400,2)] [zeros(400,2);XPM2] ones(800,1)];
    b = pinv(X)*Y; Beta3(subject,:) = b([2 4]);
    Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R23(subject) = SSeffect / SStotal;
    
    % model 4: BF + PM
    x11 = conv(X1,xBF.bf(:,1)); x11 = x11(1:400)+100;
    x12 = conv(X1,xBF.bf(:,2)); x12 = x12(1:400)+100;
    x13 = conv(X1,xBF.bf(:,3)); x13 = x13(1:400)+100;
    x2 = X1; x2(X1==1) = detrend(RT_condition1,'constant'); % parametric regressor cond1
    x21 = conv(x2,xBF.bf(:,1)); x21 = x21(1:400)+100;
    x22 = conv(x2,xBF.bf(:,2)); x22 = x22(1:400)+100;
    x23 = conv(x2,xBF.bf(:,3)); x23 = x23(1:400)+100;
    XPM1 = spm_orth([x11 x12 x13 x21 x22 x23]); % x with PM for Y1
    x2 = X1; x2(X1==1) = detrend(RT_condition2,'constant'); % parametric regressor cond2
    x21 = conv(x2,xBF.bf(:,1)); x21 = x21(1:400)+100;
    x22 = conv(x2,xBF.bf(:,2)); x22 = x22(1:400)+100;
    x23 = conv(x2,xBF.bf(:,3)); x23 = x23(1:400)+100;
    XPM2 = spm_orth([x11 x12 x13 x21 x22 x23]); % x with PM for Y2
    X = [[XPM1;zeros(400,6)] [zeros(400,6);XPM2] ones(800,1)];
    b = pinv(X)*Y; Beta4(subject,:) = b([4 10]);
    b4corrected = (sqrt(((b(4)^2)*sum(X(:,4).^2))+((b(5)^2)*sum(X(:,5).^2))+((b(6)^2)*sum(X(:,6).^2))))*sign(b(4));
    b10corrected = (sqrt(((b(10)^2)*sum(X(:,10).^2))+((b(11)^2)*sum(X(:,11).^2))+((b(12)^2)*sum(X(:,12).^2))))*sign(b(10));
    Beta4corrected(subject,:) = [b4corrected b10corrected]; Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R24(subject) = SSeffect / SStotal;
    
    % model 5: hrf mean RT
    x1 = conv(X2,xBF.bf(:,1)); x1 = x1(1:400)+100;
    X = [[x1;zeros(400,1)] [zeros(400,1);x1] ones(800,1)];
    b = pinv(X)*Y;  Beta5(subject,:) = b(1:2);
    Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R25(subject) = SSeffect / SStotal;
    
    % model 6: BF mean RT
    x11 = conv(X2,xBF.bf(:,1)); x11 = x11(1:400)+100;
    x12 = conv(X2,xBF.bf(:,2)); x12 = x12(1:400)+100;
    x13 = conv(X2,xBF.bf(:,3)); x13 = x13(1:400)+100;
    X = [[[x11 x12 x13];zeros(400,3)] [zeros(400,3);[x11 x12 x13]] ones(800,1)];
    b = pinv(X)*Y; Beta6(subject,:) = b([1 4]);
    b1corrected = (sqrt(((b(1)^2)*sum(X(:,1).^2))+((b(2)^2)*sum(X(:,2).^2))+((b(3)^2)*sum(X(:,3).^2))))*sign(b(1));
    b4corrected = (sqrt(((b(4)^2)*sum(X(:,4).^2))+((b(5)^2)*sum(X(:,5).^2))+((b(6)^2)*sum(X(:,6).^2))))*sign(b(4));
    Beta6corrected(subject,:) = [b1corrected b4corrected]; Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R26(subject) = SSeffect / SStotal;
    
    % model 7: hrf mean RT + PM
    x1 = conv(X2,xBF.bf(:,1)); x1 = x1(1:400)+100;
    PM = repmat(detrend(RT_condition1,'constant')',[1,duration])';
    x2 = X2; x2(X2==1) = PM(:); % mean center RT
    x2 = conv(x2,xBF.bf(:,1)); x2 = x2(1:400)+100;
    XPM1 = spm_orth([x1 x2]); % x with parameteric modulation for Y1
    PM = repmat(detrend(RT_condition2,'constant')',[1,duration])';
    x2 = X2; x2(X2==1) = PM(:); % mean center RT
    x2 = conv(x2,xBF.bf(:,1)); x2 = x2(1:400)+100;
    XPM2 = spm_orth([x1 x2]); % x with parameteric modulation for Y2
    X = [[XPM1;zeros(400,2)] [zeros(400,2);XPM2] ones(800,1)];
    b = pinv(X)*Y; Beta7(subject,:) = b([2 4]);
    Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R27(subject) = SSeffect / SStotal;
    
    % model 8: BF mean RT + PM
    x11 = conv(X2,xBF.bf(:,1)); x11 = x11(1:400)+100;
    x12 = conv(X2,xBF.bf(:,2)); x12 = x12(1:400)+100;
    x13 = conv(X2,xBF.bf(:,3)); x13 = x13(1:400)+100;
    PM = repmat(detrend(RT_condition1,'constant')',[1,duration])';
    x2 = X2; x2(X2==1) = PM(:); % mean center RT
    x21 = conv(x2,xBF.bf(:,1)); x21 = x21(1:400)+100;
    x22 = conv(x2,xBF.bf(:,2)); x22 = x22(1:400)+100;
    x23 = conv(x2,xBF.bf(:,3)); x23 = x23(1:400)+100;
    XPM1 = spm_orth([x11 x12 x13 x21 x22 x23]); % x with PM for Y1
    PM = repmat(detrend(RT_condition2,'constant')',[1,duration])';
    x2 = X2; x2(X2==1) = PM(:); % mean center RT
    x21 = conv(x2,xBF.bf(:,1)); x21 = x21(1:400)+100;
    x22 = conv(x2,xBF.bf(:,2)); x22 = x22(1:400)+100;
    x23 = conv(x2,xBF.bf(:,3)); x23 = x23(1:400)+100;
    XPM2 = spm_orth([x11 x12 x13 x21 x22 x23]); % x with PM for Y2
    X = [[XPM1;zeros(400,6)] [zeros(400,6);XPM2] ones(800,1)];
    b = pinv(X)*Y; Beta8(subject,:) = b([4 10]);
    b4corrected = (sqrt(((b(4)^2)*sum(X(:,4).^2))+((b(5)^2)*sum(X(:,5).^2))+((b(6)^2)*sum(X(:,6).^2))))*sign(b(4));
    b10corrected = (sqrt(((b(10)^2)*sum(X(:,10).^2))+((b(11)^2)*sum(X(:,11).^2))+((b(12)^2)*sum(X(:,12).^2))))*sign(b(10));
    Beta8corrected(subject,:) = [b4corrected b10corrected]; Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R28(subject) = SSeffect / SStotal;
    
    % model 9: hrf with durations = RT
    x1 = conv(X3,xBF.bf(:,1)); x1 = x1(1:400)+100;
    x2 = conv(X4,xBF.bf(:,1)); x2 = x2(1:400)+100;
    X = [[x1;zeros(400,1)] [zeros(400,1);x2] ones(800,1)];
    b = pinv(X)*Y; Beta9(subject,:) = b(1:2);
    Yhat = X*b;
    SSeffect = norm(Yhat-mean(Yhat)).^2;
    R29(subject) = SSeffect / SStotal;
    
end

%% 1st level fit

MeanR2 = mean([R21' R22' R23' R24' R25' R26' R27' R28' R29'],1);
RT_diff = mean(RT_diff);

%% 2nd level

% test detection
[Power(1),~,~,STATS] = ttest(Beta1(:,2)+Beta1(:,1));

[Power(2),~,~,STATS] = ttest(Beta2(:,2)+Beta2(:,1));

[Power(3),~,~,STATS] = ttest(Beta3(:,2)+Beta3(:,1));

[Power(4),~,~,STATS] = ttest(Beta4(:,2)+Beta4(:,1));

[Power(5),~,~,STATS] = ttest(Beta5(:,2)+Beta5(:,1));

[Power(6),~,~,STATS] = ttest(Beta6(:,2)+Beta6(:,1));

[Power(7),~,~,STATS] = ttest(Beta7(:,2)+Beta7(:,1));

[Power(8),~,~,STATS] = ttest(Beta8(:,2)+Beta8(:,1));

[Power(9),~,~,STATS] = ttest(Beta9(:,2)+Beta9(:,1));

[Power(10),~,~,STATS] = ttest(Beta2corrected(:,2)+Beta2corrected(:,1));

[Power(11),~,~,STATS] = ttest(Beta4corrected(:,2)+Beta4corrected(:,1));

[Power(12),~,~,STATS] = ttest(Beta6corrected(:,2)+Beta6corrected(:,1));

[Power(13),~,~,STATS] = ttest(Beta8corrected(:,2)+Beta8corrected(:,1));

Power_detection = Power;

% test the difference
[Power(1),~,~,STATS] = ttest(Beta1(:,2)-Beta1(:,1));

[Power(2),~,~,STATS] = ttest(Beta2(:,2)-Beta2(:,1));

[Power(3),~,~,STATS] = ttest(Beta3(:,2)-Beta3(:,1));

[Power(4),~,~,STATS] = ttest(Beta4(:,2)-Beta4(:,1));

[Power(5),~,~,STATS] = ttest(Beta5(:,2)-Beta5(:,1));

[Power(6),~,~,STATS] = ttest(Beta6(:,2)-Beta6(:,1));

[Power(7),~,~,STATS] = ttest(Beta7(:,2)-Beta7(:,1));

[Power(8),~,~,STATS] = ttest(Beta8(:,2)-Beta8(:,1));

[Power(9),~,~,STATS] = ttest(Beta9(:,2)-Beta9(:,1));

[Power(10),~,~,STATS] = ttest(Beta2corrected(:,2)-Beta2corrected(:,1));

[Power(11),~,~,STATS] = ttest(Beta4corrected(:,2)-Beta4corrected(:,1));

[Power(12),~,~,STATS] = ttest(Beta6corrected(:,2)-Beta6corrected(:,1));

[Power(13),~,~,STATS] = ttest(Beta8corrected(:,2)-Beta8corrected(:,1));

Power_difference = Power;

