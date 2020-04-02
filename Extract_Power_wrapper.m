function [position3,graph2] = Extract_Power_wrapper(figN,figTitle,TF,t,f)

[TF] = Extract_Power(TF,t,f);

figure(figN); clf % participants mean
graph1 = subplot(211);
% TF_O_high.clim = [-1 1]; % set colour limit of the plot
TF_Plot(mean(TF.ALL,3),TF,figTitle)
Plot_DrawPeakArea(t,f)

graph2 = subplot(212);
% plot(TF.times,mean(TF.stats.ALL.power_time,1),'LineWidth',2);
hold on
plot(TF.times,mean(TF.stats.Hits.power_time,1),'color',[217/255,83/255,25/255],'LineWidth',1)
plot(TF.times,mean(TF.stats.CR.power_time,1),'color',[126/255, 47/255, 142/255],'LineWidth',1)

% line([0,0],[-1 1],'color',[.5 .5 .5],'LineWidth',.5)
line([TF.times(1) TF.times(end)],[0 0],'Color',[.5 .5 .5],'LineStyle','--')
xlabel('Time (ms)'), ylabel('Power (dB)'), title(['Averaged power ' num2str(f(1)) '-' num2str(f(2)) 'Hz'])
xlim([-300 1200])
% legend({'All','CRs','Hits'},'Location','EastOutside')
legend({'Hits','CRs'},'Location','EastOutside')
box off
legend boxoff

% gcf(set
position3 = graph1.Position(3); 
graph2.Position(3) = position3;% to align the figures since the legend on the side distracts

% statistics: compare the power for area t by f across participants with paired t test 
[H,P,CI,~] = ttest2(TF.stats.Hits.power,TF.stats.CR.power);

disp(['H = ' num2str(H) '; p = ' num2str(P) '; CI = ' num2str(CI(1)) ' - '  num2str(CI(2))])
    
end


function [] = Plot_DrawPeakArea(t,f)
% simply draw a peak area on the tf image and write the values

x = t(1);
y = f(1);
w = t(2) - t(1);
h = f(2) - f(1);
rectangle('Position',[x y w h]) % [x y w h]

if t(2) <= 1000
    text(t(2)+10,f(2),{['T: ' num2str(t(1)) ' - ' num2str(t(2)) ' ms' ],[' F: ' num2str(f(1)) ' - ' num2str(f(2)) ' Hz']})
else
    text(t(1)-400,f(2),{['T: ' num2str(t(1)) ' - ' num2str(t(2)) ' ms' ],[' F: ' num2str(f(1)) ' - ' num2str(f(2)) ' Hz']})
end

end