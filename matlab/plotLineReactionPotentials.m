figure;
set(gca,'fontsize',16)
plot(Zgrid,analyticalReactionPotentialMap,'k','linewidth',2);
hold on
plot(Zgrid(1:4:end),bemPanelReactionPotentialMap(1:4:end),'ro','linewidth',2,'markersize',12)
plot(Zgrid(1:4:end),bemPointReactionPotentialMap(1:4:end),'bs','linewidth',2,'markersize',12)
xlabel('Field location: (0, 0, Z Angstroms)')
ylabel('Reaction potential (kcal/mol/e)')
legend('Analytical','Panel BEM','Point BEM');
plot(Zgrid,analyticalReactionPotentialMap,'k','linewidth',2);
