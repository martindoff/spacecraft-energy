L = [432 576 864 1728 3456 4320 5082 7200 7854 8640 9600 10800  12342 14400 15840 17280];
time = [0.5489  0.6366  0.8576 1.9131 5.6635  5.5380 7.8699 12.7709 20.7329 23.3681 32.0904  38.8187 54.1515 72.4321 103.4861 120.0620];
x = L(1):0.1:L(end);

loglog(L, time, 'ob', 'LineWidth', 2);

L = [4320 5082 7200 7854 8640 9600 10800  12342 14400 15840 17280];
time = [5.5380 7.8699 12.7709 20.7329 23.3681 32.0904  38.8187 54.1515 72.4321 103.4861 120.0620];
x1 = 2965:0.1:4079;
x2 = 4079:0.1:L(end);
% Fit line to data using polyfit
c = polyfit(L,time,2);
% Display evaluated equation y = m*x + b
%disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
% Evaluate fit equation using polyval
y = polyval(c,x2);
y_ = polyval(c,x1);
hold on 
loglog(x2,y, '-b', 'LineWidth', 2)
loglog(x1,y_, '--b', 'LineWidth', 2)
xlabel('Horizon N (-)', 'fontsize',15,'Interpreter','latex')
ylabel('Time to completion (s)', 'fontsize',15,'Interpreter','latex')
legend({'data', '$\mathcal{O}(N^2)$'}, 'fontsize',15,'Interpreter','latex')
grid on
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',15)


