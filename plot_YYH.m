function plot_YYH(Y,Yhat)

figure; 
subplot(1,2,1); plot(Y(1:400)); hold on; plot(Yhat(1:400),'r'); axis tight
subplot(1,2,2); plot(Y(400:800)); hold on; plot(Yhat(400:800),'r'); axis tight
