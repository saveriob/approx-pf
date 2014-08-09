figure(figno)

subplot(211)
	plot(1:(n-1), abs(u_true_nom), 'ko ', 1:(n-1), abs(u_appr_nom), 'k. ');
	hold on
	plot(1:(n-1), abs(u_true), 'ro ', 1:(n-1), abs(u_appr), 'r. ');
	title('Voltage magnitude')
	xlim([1 n-1]);
subplot(212)
	plot(1:(n-1), angle(u_true_nom), 'ko ', 1:(n-1), angle(u_appr_nom), 'k. ');
	hold on
	plot(1:(n-1), angle(u_true), 'ro ', 1:(n-1), angle(u_appr), 'r. ');
	title('Voltage angle')
	xlim([1 n-1]);
