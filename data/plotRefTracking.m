% unpack
load('ackermannData.mat');
t = ackermannData.t;
x = ackermannData.x;
xref = ackermannData.x_ref0;

% consts
r2d = 180/pi;

% plot
figure

ax1 = subplot(6,1,1);
hold on; grid on;
plot(t, xref(1,:), '--r')
plot(t, x(1,:), '-k')
title('$X$', 'interpreter', 'latex')
ylabel('[m]')
legend({'cmd', 'actual'});

ax2 = subplot(6,1,2);
hold on; grid on;
plot(t, xref(2,:), '--r')
plot(t, x(2,:), '-k')
title('$Y$', 'interpreter', 'latex')
ylabel('[m]')

ax3 = subplot(6,1,3);
hold on; grid on;
plot(t, wrapToPi(xref(3,:))*r2d, '--r')
plot(t, wrapToPi(x(3,:))*r2d, '-k')
title('$\theta$', 'interpreter', 'latex')
ylabel('[deg]')

ax4 = subplot(6,1,4);
hold on; grid on;
plot(t, xref(4,:)*r2d, '--r')
plot(t, x(4,:)*r2d, '-k')
title('$\psi$', 'interpreter', 'latex')
ylabel('[deg]')

ax5 = subplot(6,1,5);
hold on; grid on;
plot(t, xref(5,:), '--r')
plot(t, x(5,:), '-k')
title('$V$', 'interpreter', 'latex')
ylabel('[m/s]')

ax6 = subplot(6,1,6);
hold on; grid on;
plot(t, xref(6,:)*r2d, '--r')
plot(t, x(6,:)*r2d, '-k')
title('$\dot{\psi}$', 'interpreter', 'latex')
ylabel('[deg/s]')
xlabel('[sec]')

linkaxes([ax1 ax2 ax3 ax4 ax5 ax6], 'x')