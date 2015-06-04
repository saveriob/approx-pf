ps = [1 2 Inf];
nps = length(ps);

margins = zeros(n-1,nps);

figure(figno)

hold off
plot(1:(n-1), abs(u_true), 'ko ', 1:(n-1), abs(u_appr), 'k. ');
title('Voltage magnitude')
hold on

for p = 1:nps

	normp = norm(X,ps(p),'rows');
	margins(:,p) = 4 * normp * max(normp) * norm(s,holder(ps(p)))^2;

	plot(1:(n-1), abs(u_appr)-margins(:,p), num2str(p), 1:(n-1), abs(u_appr)+margins(:,p), num2str(p));
	text(n, abs(u_appr(end))-margins(end,p), num2str(holder(ps(p))));

end

xlim([1 n-1]);
ylim([0.8 1]);
