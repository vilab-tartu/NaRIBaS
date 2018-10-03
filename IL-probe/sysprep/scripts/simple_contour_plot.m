load FILE_NAME

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.9 1.9])
NAME = flipud(NAME);     %Flip matrix up to down
NAME = transpose(NAME);  %Transpose matrix
contour(NAME,'LineWidth',1.5);
set(gca,'position',[0.15 0.15 0.8 0.8],'LineWidth',0.75);
axis square
fLineWidth = 1; % width of the line of the axes
title('Contour for MOLS','FontSize',10);

ylabel('$R / \mathrm{nm}$','Interpreter','LaTeX','FontSize',10);
xlabel('$\Delta z / \mathrm{nm}$','Interpreter','LaTeX','FontSize',10);

YTick = get(gca,'YTick');
YTick = YTick/20 - 0.5;
set(gca,'YTickLabel',YTick,'FontSize',8,'LineWidth',0.75);

XTick = get(gca,'XTick');
XTick = XTick/20 - 0.5 + WALL;
set(gca,'XTickLabel',XTick,'FontSize',8,'LineWidth',0.75);
caxis([-0.012, 0.008])
;c=colorbar;
;set(c,'FontSize',10);

xpos = CENTER/2 - DIAMETER/2;
ypos = CENTER/2 - DIAMETER/2;
rectangle ('position',[xpos, ypos, DIAMETER, DIAMETER], 'curvature',[1, 1], 'LineWidth',1.5, 'LineStyle','-','FaceColor','r');
xpos = CENTER/2 - 12.5;
ypos = CENTER/2 - 12.5;
rectangle ('position',[xpos, ypos, 25, 25], 'curvature',[1, 1], 'LineWidth',1.5, 'LineStyle','--');

saveas(gcf(), 'MOLS_NAME.png', 'png');
saveas(gcf(), 'MOLS_NAME.eps', 'eps2c');
exit;
