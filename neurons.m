% exemple of neural events

duration = abs(randn(1,100))/10; % randn(1,100);
neuron1{1} = duration;
neuron1{1}(2:10:end) = 1;
neuron1{2} = zeros(1,10);
neuron1{3} = neuron1{1};

neuron2{1} = duration;
neuron2{1}(2:8:end) = 1;
neuron2{2} = zeros(1,10);
neuron2{3} = duration;
neuron2{3}(1:15:end) = 1;

neuron3{1} = duration;
neuron3{2} = zeros(1,10);
neuron3{3} = [duration abs(randn(1,20))/10];
neuron3{1}(2:10:70) = 1;
neuron3{3}(2:10:end) = 1;

figure
subplot(3,1,1); plot(cell2mat(neuron1),'LineWidth',3);
axis([-1 230 0.1 1.01]); grid on; title('Same response trial to trial','FontSize',14)
subplot(3,1,2); plot(cell2mat(neuron2),'LineWidth',3);
axis([-1 230 0.1 1.01]); grid on; title('Change in firing rate between trials','FontSize',14)
subplot(3,1,3);plot(cell2mat(neuron3),'LineWidth',3);
axis([-1 230 0.1 1.01]); grid on; title('Change in response duration between trials','FontSize',14)

