load FILE_NAME

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3.8 1.9])
NAME = flipud(NAME);     %Flip matrix up to down
NAME = transpose(NAME);  %Transpose matrix
contour(NAME,'LineWidth',1.5, 'LevelStep',0.001);
set(gca,'position',[0.15 0.15 0.8 0.8],'LineWidth',0.75);
%axis square
fLineWidth = 1; % width of the line of the axes
%title('Contour for MOLS','FontSize',10);

%ylabel('$R / \mathrm{nm}$','Interpreter','LaTeX','FontSize',10);
%xlabel('$\Delta z / \mathrm{nm}$','Interpreter','LaTeX','FontSize',10);

YTick = get(gca,'YTick');
YTick = YTick/20-1;
set(gca,'YTickLabel',YTick,'FontSize',8,'LineWidth',0.75);

XTick = get(gca,'XTick');
XTick = XTick/20;
set(gca,'XTickLabel',XTick,'FontSize',8,'LineWidth',0.75);
caxis([-0.012, 0.008]);
%c=colorbar;
%set(c,'FontSize',10);

xpos = WALL*20 - DIAMETER/2;
ypos = CENTER/2 - DIAMETER/2;
rectangle ('position',[xpos, ypos, DIAMETER, DIAMETER], 'curvature',[1, 1], 'LineWidth',1.5, 'LineStyle','-','FaceColor','r');
xpos = WALL*20 - 10;
ypos = CENTER/2 - 10;
rectangle ('position',[xpos, ypos, 20, 20], 'curvature',[1, 1], 'LineWidth',1.5, 'LineStyle','--');

saveas(gcf(), 'MOLS_NAME.png', 'png');
saveas(gcf(), 'MOLS_NAME.eps', 'eps2c');
exit;
