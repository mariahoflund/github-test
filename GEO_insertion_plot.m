

lambda = out.sync_elem.lambda.Data;
dlam = out.d_lamda.Data;
sma = out.coe.sms.Data;
d_sma = sma - A;
pos_eci = out.pos_eci.Data;
r = sqrt(pos_eci(:,1).^2 + pos_eci(:,2).^2);
dr = r - A;

%% Plot long, dr, dz

figure;
plot3(rad2deg(lambda), dr*1e-3, pos_eci(:,3)*1e-3)
hold on
slot_side = 75/2;
plot3([min(rad2deg(lambda)), max(rad2deg(lambda))], slot_side*ones(1,2), slot_side*ones(1,2), 'Color','black')
slot_side = -75/2;
plot3([min(rad2deg(lambda)), max(rad2deg(lambda))], slot_side*ones(1,2), slot_side*ones(1,2), 'Color','black')
plot3([min(rad2deg(lambda)), max(rad2deg(lambda))], -slot_side*ones(1,2), slot_side*ones(1,2), 'Color','black')
plot3([min(rad2deg(lambda)), max(rad2deg(lambda))], slot_side*ones(1,2), -slot_side*ones(1,2), 'Color','black')
xlabel('Longitude (deg)')
ylabel('\Deltar (km)')
zlabel('\Deltaz (km)')

%% Manoeuvre data

q_eci2slo = out.q_eci2slo.Data;
acc_eci = out.acc_eci.Data;
acc_slo = quatrotate(q_eci2slo, acc_eci);
mode = out.mode.Data;
dv = out.dv.Data;

figure;
subplot(3, 1, 1)
plot(mode)
yticks([1, 2, 3])
yticklabels({'Sense', 'Thrust', 'Wait'})
zoom(gca, 0.9)
subplot(3, 1, 2)
plot(dv)
zoom(gca, 0.9)
ylabel('\Deltav (m/s)')
subplot(3, 1, 3)
plot(acc_slo)
zoom(gca, 0.9)
ylabel('acc slo (m/s^2)')
legend(["r", "v", "z"])

%% Plot lambda

lambdaM = out.sync_elem.lambdaM.Data;
lambda0 = out.sync_elem.lambda0.Data;

figure;
plot(lambda)
hold on
plot(lambdaM)
plot(lambda0)
legend(["lambda", "lambdaM", "lambda0"])

%% Plot e and i

time = out.tout;
ex = out.sync_elem.ex.Data;
ey = out.sync_elem.ey.Data;
e = [ex, ey];
e_mag = vecnorm(e, 2, 2);
e_arg = angle(ex+ey*1i);

ix = out.sync_elem.ix.Data;
iy = out.sync_elem.iy.Data;
i = [ix, iy];
i_mag = vecnorm(i, 2, 2);
i_arg = angle(ix + iy*1i);

figure;
subplot(2, 1, 1)
plot(e_mag)
hold on
plot(i_mag)
legend(["e", "i"])
title('Magnitude')
subplot(2, 1, 2)
plot(rad2deg(e_arg))
hold on
plot(rad2deg(i_arg))
legend(["e", "i"])
title('Argument')


figure;
plot3(time, ex, ey)
xlabel('Time (s)')
ylabel('ex')
zlabel('ey')
hold on
plot3(time, ix, iy)
ylim([-0.01, 0.01])
zlim([-0.01, 0.01])