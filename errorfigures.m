angle_error = abs(angle(u_true) - angle(u_appr));
angle_relerror = angle_error ./ abs(angle(u_true));
ave_angle_error = mean(angle_error)/pi*180;
max_angle_error = max(angle_error)/pi*180;
ave_angle_relerror = mean(angle_relerror)*100;
max_angle_relerror = max(angle_relerror)*100;

voltage_error = abs(abs(u_true) - abs(u_appr));
voltage_relerror = voltage_error ./ (1-abs(u_true));
ave_voltage_error = mean(voltage_error);
max_voltage_error = max(voltage_error);
ave_voltage_relerror = mean(voltage_relerror)*100;
max_voltage_relerror = max(voltage_relerror)*100;

fprintf(1, 'Voltage absolute error [p.u.]:  %f (avg.), %f (max)\n', [ave_voltage_error max_voltage_error]);
fprintf(1, 'Voltage relative error [%%]:     %f (avg.), %f (max)\n', [ave_voltage_relerror max_voltage_relerror]);
fprintf(1, 'Angle absolute error [p.u.]:    %f (avg.), %f (max)\n', [ave_angle_error max_angle_error]);
fprintf(1, 'Angle relative error [%%]:       %f (avg.), %f (max)\n', [ave_angle_relerror max_angle_relerror]);
