clear all
%close all
clc

addpath('matpower4.1');
define_constants;

% GRID MODEL
mpc = loadcase('case_ieee123_PQ');
%mpc = loadcase('case118');
%mpc = loadcase('case30');

PCCindex = find(mpc.bus(:,BUS_TYPE)==3);
n = length(mpc.bus(:,BUS_TYPE));
PQnodes = setdiff(1:n,PCCindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

nbr = size(mpc.branch,1);
nbu = size(mpc.bus,1);
L = zeros(nbu,nbu);

for br = 1:nbr
	br_F_BUS = mpc.branch(br,F_BUS);
	br_T_BUS = mpc.branch(br,T_BUS);
	br_BR_R = mpc.branch(br,BR_R);
	br_BR_X = mpc.branch(br,BR_X);
	br_Y = 1 / (br_BR_R + 1j * br_BR_X);

	L(br_F_BUS, br_T_BUS) = br_Y;
	L(br_T_BUS, br_F_BUS) = br_Y;
	L(br_F_BUS, br_F_BUS) = L(br_F_BUS, br_F_BUS) - br_Y;
	L(br_T_BUS, br_T_BUS) = L(br_T_BUS, br_T_BUS) - br_Y;
end

X = inv(L(PQnodes,PQnodes));

e0 = zeros(n,1);
e0(PCCindex) = 1;
X2 = inv([L,e0;e0',0]); X2 = X2(1:n,1:n);
X2 = X2(PQnodes,PQnodes);

%%%%%%%

s = mpc.bus(PQnodes,PD) + mpc.bus(PQnodes,GS) + 1j * (mpc.bus(PQnodes,QD) - mpc.bus(PQnodes,BS));

rho_2_star = max(norm(X,2,'cols'));
fprintf(1,'X 2-star: %f\n', rho_2_star);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2_star * norm(s,2));

rho_2 = norm(X,2);
fprintf(1,'X 2-norm: %f\n', rho_2);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2 * norm(s,2));

rho_1_star = max(max(abs(X)));
fprintf(1,'X 1-star: %f\n', rho_1_star);
fprintf(1,'s 1-norm: %f\n', norm(s,1));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_1_star * norm(s,1));

rho_Inf_star = norm(X, Inf);
fprintf(1,'X inf-star: %f\n', rho_Inf_star);
fprintf(1,'s inf-norm: %f\n', norm(s,Inf));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_Inf_star * norm(s,Inf));

%%%%%%%%%

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);
u_tild = abs(u_true - u_appr);

f = X\(u_true - 1 - X * conj(s));

%%%%%%%%%%

ps = [1 2 Inf];
nps = length(ps);

margins = zeros(n-1,nps);

figure(1)

hold off
plot(1:(n-1), results.bus(PQnodes,VM), 'ko ', 1:(n-1), abs(1 + X * conj(s)), 'k. ');
title('Voltage magnitude')
hold on

D = eye(n-1);
D(112) = 100;

for p = 1:nps

	normp = norm(X,ps(p),'rows');
	margins(:,p) = 4 * normp * max(normp) * norm(s,holder(ps(p)))^2;

	plot(1:(n-1), abs(u_appr)-margins(:,p), num2str(p), 1:(n-1), abs(u_appr)+margins(:,p), num2str(p));
	text(n, abs(u_appr(end))-margins(end,p), num2str(holder(ps(p))));

end

%%%%%%%%%%

% DC power flow comparison

figure(2)

plot(1:(n-1), results.bus(PQnodes,VA), 'ko ', 1:(n-1), angle(1 + X * conj(s))/pi*180, 'k. ');
title('Voltage angle')

hold on

% linear angles
plot(1:(n-1), imag(X * conj(s))/pi*180, 'kx ');

% DC power flow
plot(1:(n-1), (imag(X) * real(s)) /pi*180, 'k+ ');

hold off



%%%%%%%%%%%

if false

figure(5)
scatter(real(u_appr),imag(u_appr),'r')
hold on
title('APPR (red) - TRUE (black) - DC (blue)')
scatter(real(u_true),imag(u_true),'k')

theta_dc = imag(X) * real(s);
u_dc = exp(j*theta_dc);

scatter(real(u_dc),imag(u_dc),'b')

xlim([0.97 1.005])
ylim([-0.05 0.001])
scatter(1,0,'k')
axis equal
hold off

end




%%%%%%%%%%%%

% Export magnitudes

if export_magnitues=true

	fname = 'data_voltagem-1norm.data';

	myfile=fopen(fname,"w");
	fdisp(myfile,'bus appr real min max');
	fclose(myfile);

	data_voltagem = [(1:n-1)'					...
			 abs(u_appr)						...
			 results.bus(PQnodes,VM)			...
			 abs(u_appr)-margins(:,3)			...
			 abs(u_appr)+margins(:,3)];

	save('-append', '-ascii', fname, 'data_voltagem');

end

% Export angles

if export_angles=false

	fname = 'data_voltagea.data';

	myfile=fopen(fname,"w");
	fdisp(myfile,'bus real approx linear dc');
	fclose(myfile);

	data_voltagea = [ ...
			(1:n-1)' ...
			results.bus(PQnodes,VA) ...
			angle(1 + X * conj(s))/pi*180 ...
			imag(X * conj(s))/pi*180 ...
			(imag(X) * real(s))/pi*180 ];

	save('-append', '-ascii', fname, 'data_voltagea');

end
