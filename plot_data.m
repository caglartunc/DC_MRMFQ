% Plotting script that plots the buffer distribution

B_limit = B;
B_gran = 5;

subplot(1,2,1)
set(gca,'FontSize',18) 
area(0:delta:B, fs_idle')
title(['(',title_text(config.params.subplot),')'],'Interpreter','latex')
xlim([0,B_limit])
ylim([0 ceil(max(sum(fs_idle,1)) * 10) / 10])
xticks(0:B_gran:B_limit);
xlabel('Buffer level $x$','Interpreter','latex')
ylabel('$\hat{f}_i(x)$','Interpreter','latex')
set(gca,'FontSize',18)
grid on
hold on
plot([Ta,Ta],[0,ceil(max(sum(fs_idle,1)) * 10) / 10],'g-.','LineWidth',2)
plot([Td,Td],[0,ceil(max(sum(fs_idle,1)) * 10) / 10],'k--','LineWidth',1)

subplot(1,2,2)
set(gca,'FontSize',18)
area(0:delta:B, Fs_idle')
title(['(',title_text(config.params.subplot+4),')'],'Interpreter','latex')
xlim([0,B_limit])
ylim([0 1])
xticks(0:B_gran:B_limit);
xlabel('Buffer level $x$','Interpreter','latex')
ylabel('$\hat{F}_i(x)$','Interpreter','latex')
set(gca,'FontSize',18) 
grid on
hold on
plot([Ta,Ta],[0,1],'g-.','LineWidth',2)
plot([Td,Td],[0,1],'k--','LineWidth',1)
if(config.params.subplot~=inf)
    legend(config.params.legend_text,'Interpreter','latex','Orientation','horizontal')
end