%% This is a simple simulation to show what happens to the estimated BOLD 
%% amplitude when the data are not well modelled, and how derivatives can 
%% remedee this and then reestimate the amplitude


%% simulate data to look at the effect of 'shape'

xBF.dt = 0.5;
xBF.name = 'hrf (with time and dispersion derivative)';
xBF.length = 32; 
xBF.order = 1;
xBF.T = 30;
xBF = spm_get_bf(xBF);
xBF.UNITS = 'secs';

onsets = [1:60:600];
Y  = zeros(700,1);
XX = zeros(700,1);
for i=1:10
    XX(onsets(i)) = 1; % design matrix with duration 0
    Y(onsets(i)) = 1;
end
Y1 = conv(Y,xBF.bf*[1 0 0]'); Y1 = Y1(1:600)./max(Y1(:))+100;
Y2 = conv(Y,xBF.bf*[1 -2 0]'); Y2 = Y2(1:600)./max(Y2(:))+100;
Y3 = conv([[0 0]';Y(1:end-2)],xBF.bf*[1 0 -1]'); Y3 = Y3(1:600)./max(Y3(:))+100;

% the data Y1 are the standard hrf, Yw2 is delayed and Y3 has a different
% shape - data are normalized by the max so they peak at value 101 ; the
% baseline will be however a little different because undershoots differ,
% in particular for Y2
plot(Y1,'LineWidth',3); hold on; plot(Y3,'g','LineWidth',3);
plot(Y2,'r','LineWidth',3); grid on; title('Simulated data','FontSize',14);

%% data fit
% now if we apply the same model to each set we can see that model
% fits Y1, for Y2 it misses the max because of the temporal delay, and for
% Y3 overshoots to fit as much as possible the larger hrf

x1 = conv(XX,xBF.bf(:,1)); x1 = x1(1:600);
X = [x1 ones(600,1)];
b = pinv(X)*Y1; Yhat = X*b; 

figure('Name','hrf model');
subplot(3,1,1); plot(Y1,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
title(['estimated b hrf =' num2str(b(1))]);

b = pinv(X)*Y2; Yhat = X*b; 
subplot(3,1,2); plot(Y2,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
title(['estimated b hrf =' num2str(b(1))]);

b = pinv(X)*Y3; Yhat = X*b; 
subplot(3,1,3); plot(Y3,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
title(['estimated b hrf =' num2str(b(1))]);

% use 1st derivative
x1 = conv(XX,xBF.bf(:,1)); x1 = x1(1:600);
x2 = conv(XX,xBF.bf(:,2)); x2 = x2(1:600);
X = [spm_orth([x1 x2]) ones(600,1)];
b = pinv(X)*Y1; Yhat = X*b; 

figure('Name','hrf + 1st deriv model');
subplot(3,1,1); plot(Y1,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
title(['estimated b hrf =' num2str(b(1))]);

b = pinv(X)*Y2; Yhat = X*b; 
subplot(3,1,2); plot(Y2,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
B = sqrt(((b(1).^2).*sum(X(:,1).^2)) + ((b(2).^2).*sum(X(:,2).^2))); 
title(['estimated b hrf =' num2str(b(1)) ' corrected =' num2str(B)]);

b = pinv(X)*Y3; Yhat = X*b; 
subplot(3,1,3); plot(Y3,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
B = sqrt(((b(1).^2).*sum(X(:,1).^2)) + ((b(2).^2).*sum(X(:,2).^2))); 
title(['estimated b hrf =' num2str(b(1)) ' corrected =' num2str(B)]);

% use 2nd derivatives
x1 = conv(XX,xBF.bf(:,1)); x1 = x1(1:600);
x2 = conv(XX,xBF.bf(:,2)); x2 = x2(1:600);
X = [spm_orth([x1 x2]) ones(600,1)];
b = pinv(X)*Y1; Yhat = X*b; 

figure('Name','hrf + 2nd deriv model');
subplot(3,1,1); plot(Y1,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
title(['estimated b hrf =' num2str(b(1))]);

b = pinv(X)*Y2; Yhat = X*b; 
subplot(3,1,2); plot(Y2,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
B = sqrt(((b(1).^2).*sum(X(:,1).^2)) + ((b(2).^2).*sum(X(:,2).^2))); 
title(['estimated b hrf =' num2str(b(1)) ' corrected =' num2str(B)]);

b = pinv(X)*Y3; Yhat = X*b; 
subplot(3,1,3); plot(Y3,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
B = sqrt(((b(1).^2).*sum(X(:,1).^2)) + ((b(2).^2).*sum(X(:,2).^2))); 
title(['estimated b hrf =' num2str(b(1)) ' corrected =' num2str(B)]);

% use all derivatives
x1 = conv(XX,xBF.bf(:,1)); x1 = x1(1:600);
x2 = conv(XX,xBF.bf(:,2)); x2 = x2(1:600);
x3 = conv(XX,xBF.bf(:,3)); x3 = x3(1:600);
X = [spm_orth([x1 x2 x3]) ones(600,1)];
b = pinv(X)*Y1; Yhat = X*b; 

figure('Name','hrf + all deriv model');
subplot(3,1,1); plot(Y1,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
title(['estimated b hrf =' num2str(b(1))]);

b = pinv(X)*Y2; Yhat = X*b; 
subplot(3,1,2); plot(Y2,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
B = sqrt(((b(1).^2).*sum(X(:,1).^2)) + ((b(2).^2).*sum(X(:,2).^2)) + ((b(3).^2).*sum(X(:,3).^2))); 
title(['estimated b hrf =' num2str(b(1)) ' corrected =' num2str(B)]);

b = pinv(X)*Y3; Yhat = X*b; 
subplot(3,1,3); plot(Y3,'LineWidth',3); 
hold on; plot(Yhat,'r','LineWidth',3); grid on;
B = sqrt(((b(1).^2).*sum(X(:,1).^2)) + ((b(2).^2).*sum(X(:,2).^2)) + ((b(3).^2).*sum(X(:,3).^2))); 
title(['estimated b hrf =' num2str(b(1)) ' corrected =' num2str(B)]);


