function RT = OHBM2015_generateRT

%% RT distribution -- we take here RT for decisions times
N = 500;
pop_condition1 = NaN(1,N);
pop_condition2 = NaN(1,N);
 for event = 1:N % draw 100 random data from a gamma distribution
    if mod(event,2) == 0
    pop_condition1(event) = gamrnd(5,10)+10;
    else
    pop_condition2(event) = gamrnd(5,10)+20; % cond1/cond2 differ 
    end
end
pop_condition1(isnan(pop_condition1))=[];
pop_condition2(isnan(pop_condition2))=[];

% rescale 
scaling = max([pop_condition1 pop_condition2])/6;
pop_condition1 = pop_condition1./scaling;
pop_condition2 = pop_condition2./scaling;
pop_condition1(pop_condition1<1) = [];
pop_condition2(pop_condition2<1) = [];

% figure; 
% set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', ...
%     'outerposition',[0 0 1 1])
% [h1,x1]=hist(pop_condition1,20);
% [h2,x2]=hist(pop_condition2,20);
% bar(x1, h1); hold on; bar(x2, h2, 'r');
% h=findobj(gca, 'Type', 'patch');
% set (h, 'FaceAlpha', 0.7); grid on
% title('populations to draw RT from','FontSize',14);

%% generate 20 subjects

subjects = 0;
add_plot = 0;
RT = NaN(2,10,20); % store samples in a 2 conditions * 10 values * 20 subjects matrix
while subjects<20
    
    % keep only 10 data points
    RT_condition1 = pop_condition1(randi(length(pop_condition1),10,1))-1;
    RT_condition2 = pop_condition2(randi(length(pop_condition2),10,1))-1;

    % check we have nice differences between conditions
    [H,P,CI,STATS] = ttest(RT_condition2-RT_condition1);
    
    if P<0.05
%         subplot(1,2,2);boxplot([RT_condition1' RT_condition2'])
%         D = mean(RT_condition1)-mean(RT_condition2);
%         title(sprintf('Exemple of a sample \n difference=%g CI [%g %g] \n t(%g)=%g p=%g' ...
%             ,D,CI(1),CI(2),STATS.df,STATS.tstat,P),'FontSize',14);
%        grid on; add_plot = 1;
        subjects = subjects+1;
        RT(1,:,subjects) = RT_condition1;
        RT(2,:,subjects) = RT_condition2;
    end
end
