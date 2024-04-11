%% fig 1.1 violin plots  
clear;clc
load all_all_RT.mat;
addpath('C:\toolboxes\gramm-master\gramm-master')
% RT.congruency=mat2cell(all_all_RT(:,5),size(all_all_RT,1),1)
RT.congruency={all_all_RT(:,5)}
RT.alertness={all_all_RT(:,7)}
RT.reactiontime=all_all_RT(:,2)
RT.ID={all_all_RT(:,6)}
%%
RT.alertness1={}
for camb= 1:size(RT.alertness{1,1},1)
    if RT.alertness{1,1}(camb)==1
       RT.alertness1(camb)= {'alert'}
     
    end
    if RT.alertness{1,1}(camb)==2
       RT.alertness1(camb)= {'drowsy'}
     
    end
end;RT.alertness1=RT.alertness1'
%%
RT.congruency1={}
for camb= 1:size(RT.congruency{1,1},1)
    if RT.congruency{1,1}(camb)==1
       RT.congruency1(camb)= {'congruent'}
     
    end
    if RT.congruency{1,1}(camb)==2
       RT.congruency1(camb)= {'incongruent'}
     
    end
end;RT.congruency1=RT.congruency1'

% Create a gramm object, provide x (year of production) and y (fuel economy) data,
% color grouping data (number of cylinders) and select a subset of the data
% g=gramm('x',RT.alertness1,'y',RT.reactiontime,'color',RT.congruency1);
%%% 
g=gramm('x',RT.congruency1,'y',RT.reactiontime,'color',RT.alertness1);

% %Violin plots
g(1,1).stat_violin('fill','transparent');
g(1,1).set_title('stat_violin()');
% 
%These functions can be called on arrays of gramm objects
g.set_names('x','Congruency','y','Reaction time','color','Alertness');
g.set_title('Reaction time');
% g(1,1).stat_violin('normalization','area','dodge',0,'fill','edge');
g(1,1).stat_boxplot('width',0.15);
% g(1,1).set_title('with stat_boxplot()');
% g(1,1).set_color_options('map','brewer_dark');
gf = copy(g);
g.axe_property('YLim',[0 2100]);
figure('Position',[100 100 800 550]);
g.draw();

% figure(1)
% hold on 
% set(gca,'ylim',[0 2600])
% hold on;
% set(gca,'xlim',[0 q10])

disp('Om')
saveas(gcf,['RT_cong_incong_alert_drowsy_violin.fig']) 
saveas(gcf,['RT_cong_incong_alert_drowsy_violin.png']) 
save2pdf('RT_cong_incong_alert_drowsy_violin',gcf,600);