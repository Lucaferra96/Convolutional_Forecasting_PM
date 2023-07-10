function []=createScatter(x,y,labelx,labely,corr_best, nrmse_best)
%create scatter

if (length(x) == length(y))
    figure;
    plot(x,y,'r*')
    grid on
    hold on
    
    plot([0 max(max([x y]))],[0 max(max([x y]))],'b--');
    ylim([0, max(y, [], 'all')]);

%     plot([0 max(max([x y]))/2],[0 max(max([x y]))],'b--');
%     plot([0 max(max([x y]))],[0 max(max([x y]))/2],'b--');

    s1=xlabel(labelx);
    s2=ylabel(labely);
    s3=text(0.02*max(max([x y])),0.95*max(max([x y])),strcat('corr= ',num2str(corr_best, '%10.3f')));
    s4=text(0.02*max(max([x y])),0.85*max(max([x y])),strcat('nrmse= ',num2str(nrmse_best, '%10.3f')));
    s=[s1 s2 s3 s4];
    set(s, 'fontsize', 25);
    set(gca, 'fontsize', 25);
end
end