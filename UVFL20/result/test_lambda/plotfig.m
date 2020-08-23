load('ALS_result');
load('Brid_result');
load('AMM_result');

nr = 1500;   nc = 1500;

SR = 0.2;   sigma = 0.1;

rstar = 10;

nuAMM  = [ 3    3.5    4      4.2    4.5    4.7     5     5.2    5.5    5.8     6     6.2    6.5    6.7    7     7.2    7.5   7.8   8 ]*10; 
 
nuMAPM = [0.5   0.6    0.7    0.8    0.9     1     1.1    1.2    1.3    1.4    1.5    1.6    1.7    1.8    1.9    2     2.1   2.2  2.3 ]*10 ;    

nuALS  = [0.01  0.05   0.1    0.2    0.3    0.4    0.5     0.6   0.7    0.8     0.9     1     1.2   1.5     1.8    2    2.4    2.8    3]; 

%% ***************************************************************

subplot(1,3,1);
[AX,H1,H2] = plotyy(nuAMM,AMM_averelerr,nuAMM,AMM_averank);
title('Algorithm 1')
legend([H1,H2],'RE','rank');
xlabel('$c_{\lambda}$','interpreter','latex','fontsize',16)
set(get(AX(1),'Ylabel'),'String','RE (%)','fontsize',16); 
set(get(AX(2),'Ylabel'),'String','rank','fontsize',16);
set(H1,'linestyle','--','marker','o','color','r','Linewidth',1.5);
set(H2,'linestyle','-','marker','s','color','b','Linewidth',1.5);
set(AX,'XColor','k','YColor','k','fontsize',16);
set(AX,'Ycolor','k','FontSize',16)

set(AX,'Xlim',[30,80])
set(AX(1),'ylim',[0.03,0.71],'ytick',[0.03:0.1:0.71])
set(AX(2),'ylim',[5,20],'ytick',[5:3:20])
grid on
hold on;

 
subplot(1,3,2);
[AX,H1,H2] = plotyy(nuMAPM,Brid_averelerr,nuMAPM,Brid_averank);
title('Algorithm 3')
legend([H1,H2],'RE','rank');
xlabel('$c_{\lambda}$','interpreter','latex','fontsize',16)
set(get(AX(1),'Ylabel'),'String','RE (%)','fontsize',16); 
set(get(AX(2),'Ylabel'),'String','rank','fontsize',16);
set(H1,'linestyle','--','marker','o','color','r','Linewidth',1.5);
set(H2,'linestyle','-','marker','s','color','b','Linewidth',1.5);
set(AX,'XColor','k','YColor','k','fontsize',16);
set(AX,'Ycolor','k','FontSize',16)

set(AX,'Xlim',[5,23])
set(AX(1),'ylim',[0.03,0.78],'ytick',[0.03:0.1:0.78])
set(AX(2),'ylim',[6,20],'ytick',[6:2:20])
grid on
hold on;



subplot(1,3,3);
[AX,H1,H2] = plotyy(nuALS,ALS_averelerr,nuALS,ALS_averank);
title('ALS')
legend([H1,H2],'RE','rank');
xlabel('$c_{\lambda}$','interpreter','latex','fontsize',16)
set(get(AX(1),'Ylabel'),'String','RE (%)','fontsize',16); 
set(get(AX(2),'Ylabel'),'String','rank','fontsize',16);
set(H1,'linestyle','--','marker','o','color','r','Linewidth',1.5);
set(H2,'linestyle','-','marker','s','color','b','Linewidth',1.5);
set(AX,'XColor','k','YColor','k','fontsize',16);
set(AX,'Ycolor','k','FontSize',16)

set(AX,'Xlim',[0.01,3])
set(AX(1),'ylim',[0.6,0.88],'ytick',[0.6:0.05:0.88])
set(AX(2),'ylim',[10,150],'ytick',[10:20:150])
grid on
hold on;










