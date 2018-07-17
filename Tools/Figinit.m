function[hf,hp] = Figinit(hf,target)
  % Set figure
    set(hf,'units','normalized','outerposition',[0 0 1 1]);
    set(hf,'menubar','none')
    set(hf,'color','w');
    hp = patch(0,0,'g','Edgecolor','g');
    hold on;
    set(gca,'xlim',target*[-1.2 1.2],'ylim',target*[-1.2 1.2]);
    set(gca,'TickLength',[0 0]) 
    axis square;
    grid on;
    box on;
end