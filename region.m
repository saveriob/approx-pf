clear all
close all
clc

addpath('matpower4.1');
define_constants;

% GRID MODEL
mpc = loadcase('case_2nodes');

n = length(mpc.bus(:,BUS_TYPE));
PCCindex = 1;
nodex = 2;
nodey = 3;
PQnodes = setdiff(1:n,PCCindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

px_min = 0;
px_max = 1.2;
py_min = 0;
py_max = 0.8;

p_step = 1e-2;

pxs = px_min:p_step:px_max;
pys = py_min:p_step:py_max;

nx = length(pxs);
ny = length(pys);

if load_feasible = true

	load feasible.mat

else

	thx = 0.1;
	thy = 1;

	feasible = zeros(nx,ny);

	for xi = 1:nx
	for yi = 1:ny

		mpc.bus(nodex,PD) = cos(thx) * pxs(xi);
		mpc.bus(nodex,QD) = sin(thx) * pxs(xi);
		mpc.bus(nodey,PD) = cos(thy) * pys(yi);
		mpc.bus(nodey,QD) = sin(thy) * pys(yi);
		results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));
		feasible(xi,yi) = results.success;

	end
	end

end

contour(pxs, pys, feasible', [0.5 0.5], 'LineWidth', 3)
axis equal
hold on

%%%%%%%%%%
% build X matrix

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

%%%%%%%

%
%	|s|_p < V_0^2 / (4 * |X|_q^*)
%

%%%%%%% p = 1 , q = Inf

X_inf_star = max(max(abs(X)));
b_1 = 1/(4*X_inf_star);

line([0 b_1], [b_1 0], 'LineWidth', 3);

%%%%%%% p = Inf , q = 1

X_1_star = norm(X, Inf);
b_inf = 1/(4*X_1_star);

line([0 b_inf b_inf], [b_inf b_inf 0], 'LineWidth', 3);

%%%%%%% p = 2 , q = 2

X_2_star = max(norm(X,2,'cols'));
b_2 = 1/(4*X_2_star);

t = 0:0.01:pi/2;

xt = cos(t) * b_2;
yt = sin(t) * b_2;
plot(xt, yt, 'k', 'LineWidth', 3);

%%%%%%%% weighted norm

wnorm = zeros(nx,ny);

Xa = abs(X);

D = diag(diag((X)));

%D = diag([0.5 1]);

for xi = 1:nx
for yi = 1:ny

	s = [pxs(xi); pys(yi)];
	wnorm(xi,yi) = norm(D*s,1) < 1/( 4 * max(max(abs(X*inv(D)))) );
	% wnorm(xi,yi) = norm(D*s,2) < 1/( 4 * max(norm(X*inv(D),2,'cols')) );
	% wnorm(xi,yi) = norm(D*s,Inf) < 1/( 4 * norm(X*inv(D), Inf) );

end
end


contour(pxs, pys, wnorm', [0.5 0.5], 'LineWidth', 3)
axis equal




return





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
