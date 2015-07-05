
leftstart = 0.06;
topstart = 0.54;
rightstart = 0.55;
botstart = 0.09;
imwidth = 0.30;
imheight = 0.45;

figure;
subplot(2,2,1);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid,contourColors, ...
         'linestyle', 'none');
axis equal; caxis([minPot maxPot]); colorbar;
set(gca,'fontsize',16,'position',[leftstart topstart imwidth imheight]);
title('(a)');

subplot(2,2,2);
contourf(Ygrid,Zgrid,rPotPointBemGrid,contourColors, ...
         'linestyle', 'none');
axis equal; caxis([minPot maxPot]); colorbar;
set(gca,'fontsize',16,'position',[rightstart topstart imwidth imheight]);
title('(b)');

subplot(2,2,3);
contourf(Ygrid,Zgrid,rPotSimpleGrid,contourColors, ...
         'linestyle', 'none');
axis equal; caxis([minPot maxPot]); colorbar;
set(gca,'fontsize',16,'position',[leftstart botstart imwidth imheight]);
title('(c)');

subplot(2,2,4);
contourf(Ygrid,Zgrid,rPotPanelBemGrid,contourColors, ...
         'linestyle', 'none');
axis equal; caxis([minPot maxPot]); colorbar;
set(gca,'fontsize',16,'position',[rightstart botstart imwidth imheight]);
title('(d)');

%%%%%%%%%%%%%%%%%%%%
errColors = [analyticalReactionPotentialMap - ...
	     bemPointReactionPotentialMap;
	     analyticalReactionPotentialMap - ...
	     bemSimpleReactionPotentialMap;
	     analyticalReactionPotentialMap- ...
	     bemPanelReactionPotentialMap];

minErrPot = min(errColors); maxErrPot = max(errColors);
errContourColors = minErrPot:(maxErrPot-minErrPot)/20:maxErrPot;

figure;
axisLimits = [-4 4 2 7];
deltaRplot = 0.25;
subplot(1,3,1);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid-rPotPointBemGrid,errContourColors,'linestyle','none');
caxis([minErrPot maxErrPot]);colorbar('southoutside');
set(gca,'fontsize',16);
axis equal;
%axis([-(R+deltaRplot) R+deltaRplot -(R+deltaRplot) R+deltaRplot]); 
axis(axisLimits);
title('(a)');

subplot(1,3,2);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid-rPotSimpleGrid,errContourColors, ...
         'linestyle', 'none');
 caxis([minErrPot maxErrPot]);colorbar('southoutside');
set(gca,'fontsize',16);
axis equal;
%axis([-(R+deltaRplot) R+deltaRplot -(R+deltaRplot) R+deltaRplot]); 
axis(axisLimits);
title('(b)');

subplot(1,3,3);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid-rPotPanelBemGrid,errContourColors, ...
         'linestyle', 'none');
caxis([minErrPot maxErrPot]);colorbar('southoutside');
set(gca,'fontsize',16);
axis equal; 
%axis([-(R+deltaRplot) R+deltaRplot -(R+deltaRplot) R+deltaRplot]); 
axis(axisLimits);
title('(c)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid,contourColors, ...
         'linestyle', 'none');
axis equal; caxis([minPot maxPot]); colorbar;
set(gca,'fontsize',16,'position',[leftstart topstart imwidth imheight]);
title('(a)');

subplot(2,2,2);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid-rPotPointBemGrid,errContourColors, ...
         'linestyle', 'none');
axis equal; axis(axisLimits);
caxis([minErrPot maxErrPot]); colorbar;
set(gca,'fontsize',16,'position',[rightstart topstart imwidth imheight]);
title('(b)');

subplot(2,2,3);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid-rPotSimpleGrid,errContourColors, ...
         'linestyle', 'none');
axis equal; axis(axisLimits)
caxis([minErrPot maxErrPot]); colorbar;
set(gca,'fontsize',16,'position',[leftstart botstart imwidth imheight]);
title('(c)');

subplot(2,2,4);
contourf(Ygrid,Zgrid,rPotAnalyticalGrid-rPotPanelBemGrid,errContourColors, ...
         'linestyle', 'none');
axis equal; axis(axisLimits)
caxis([minErrPot maxErrPot]); colorbar;
set(gca,'fontsize',16,'position',[rightstart botstart imwidth imheight]);
title('(d)');

print -depsc2 analytical-vs-bem-errors.eps