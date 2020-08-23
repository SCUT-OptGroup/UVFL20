load('scheme1_result');
load('scheme2_result');

%% scheme 1 :  nuAMM  = 60;    nuMAPM = 15;   nuALS = 1; 

%% scheme 2 :  nuAMM  = 66.8;  nuMAPM = 15;   nuALS = 1;

nr = 1000;  nc = 1000;

rstar = 5;   noiseratio = 0.1;

SR_list = [0.04   0.06   0.08    0.10     0.12    0.14   0.16    0.18   0.2];

%% ***************************************************************************************************
subplot(1,2,1);
h=plot(SR_list,AMM_averelerr1,'rs-', SR_list,Brid_averelerr1,'b*-.',SR_list,ALS_averelerr1,'md-',SR_list,Toh_averelerr1,'g+:');
set(h,'LineWidth',1.5) 
xlabel('Sampling Ratio ');   ylabel('RE');
set(get(gca,'XLabel'),'FontSize',14);
set(get(gca,'YLabel'),'FontSize',14);
legend('Algorithm 1','Algorithm 3','ALS','ADMM');
title('Relative Error for (r,n)=(5,1000) under Scheme 1')
grid on
hold on;


subplot(1,2,2);
h=plot(SR_list,AMM_averelerr2,'rs-', SR_list,Brid_averelerr2,'b*-.',SR_list,ALS_averelerr2,'md-',SR_list,Toh_averelerr2,'g+:');
set(h,'LineWidth',1.5) 
xlabel('Sampling Ratio ');   ylabel('RE');
set(get(gca,'XLabel'),'FontSize',14);
set(get(gca,'YLabel'),'FontSize',14);
legend('Algorithm 1','Algorithm 3','ALS','ADMM');
title('Relative Error for (r,n)=(5,1000) under Scheme 2')
grid on
hold on;

















