clear all
close all
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

tot = sum(mpc.bus(PQnodes,PD)) + sum(mpc.bus(PQnodes,GS)) + 1j * (sum(mpc.bus(PQnodes,QD)) - sum(mpc.bus(PQnodes,BS)));
mpc.bus(PQnodes,PD) = 0;
mpc.bus(PQnodes,GS) = 0;
mpc.bus(PQnodes,QD) = 0;
mpc.bus(PQnodes,BS) = 0;

mpc.bus(114,PD) = real(tot);
mpc.bus(114,QD) = imag(tot);

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));
%results = runpf(mpc);

%v = results.bus(:,VM);
%losses = sum(real(get_losses(results)));

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

%rho_2_star = sqrt(0.030469);
rho_2_star = sqrt(0.039461);
fprintf(1,'X 2-star: %f\n', rho_2_star);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2_star * norm(s,2));

rho_2_col = max(norm(X,2,'cols'));
fprintf(1,'X max col norm: %f\n', rho_2_col);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2_col * norm(s,2));

rho_2 = norm(X,2);
fprintf(1,'X 2-norm: %f\n', rho_2);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_2 * norm(s,2));

rho_1_max = max(max(abs(X)));
fprintf(1,'X max-norm: %f\n', rho_1_max);
fprintf(1,'s 1-norm: %f\n', norm(s,1));
fprintf(1,'4 * ||X|| * ||s|| = %f\n\n', rho_1_max * norm(s,1));

%%%%%%

figure(1)
plot(1:(n-1), results.bus(PQnodes,VM), 'ko ', 1:(n-1), abs(1 + X * conj(s)), 'k. ');
title('Voltage magnitude')

figure(2)
plot(1:(n-1), results.bus(PQnodes,VA), 'ko ', 1:(n-1), angle(1 + X * conj(s))/pi*180, 'k. ');
title('Voltage angle')


%%%%%%%%%

margin_norm1 = 4 * diag(abs(X)) * rho_1_max * norm(s,1)^2;
margin_rho2col = 4 * norm(X,2,'cols')' * rho_2_col * norm(s,2)^2;
margin_tightest = 4 * norm(X,2,'cols')' * rho_2_star * norm(s,2)^2;
margin_uniform = 4 * rho_2_col * rho_2_star * norm(s,2)^2;

figure(1);
hold on
%plot(1:(n-1), results.bus(PQnodes,VM)-margin_tightest, 'k-', 1:(n-1), results.bus(PQnodes,VM)+margin_tightest, 'k');
plot(1:(n-1), results.bus(PQnodes,VM)-margin_norm1,    'b-', 1:(n-1), results.bus(PQnodes,VM)+margin_norm1,    'b');
%plot(1:(n-1), results.bus(PQnodes,VM)-margin_uniform, 'r-', 1:(n-1), results.bus(PQnodes,VM)+margin_uniform, 'r');

figure(2);
hold on
plot(1:(n-1), results.bus(PQnodes,VA)-margin_tightest/pi*180, 'k', 1:(n-1), results.bus(PQnodes,VA)+margin_tightest/pi*180, 'k');
%plot(1:(n-1), results.bus(PQnodes,VA)-margin_norm1/pi*180, 'b', 1:(n-1), results.bus(PQnodes,VA)+margin_norm1/pi*180, 'b');
%plot(1:(n-1), results.bus(PQnodes,VA)-margin_uniform/pi*180, 'r', 1:(n-1), results.bus(PQnodes,VA)+margin_uniform/pi*180, 'r');


u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);
u_tild = abs(u_true - u_appr);

f = X\(u_true - 1 - X * conj(s));

%figure(3)
%plot(1:(n-1), u_tild, 'k* ', 1:(n-1), margin_rho2col, 'k-');

%%%%%%%%%%%%

% Export

fname = 'data_voltagem.data';

myfile=fopen(fname,"w");
fdisp(myfile,'bus appr real min max');
fclose(myfile);

data_voltagem = [(1:n-1)' abs(1 + X * conj(s)) results.bus(PQnodes,VM) results.bus(PQnodes,VM)-margin_norm1 results.bus(PQnodes,VM)+margin_norm1];

save('-append', '-ascii', fname, 'data_voltagem');




fname = 'data_voltagea.data';

myfile=fopen(fname,"w");
fdisp(myfile,'bus appr real min max');
fclose(myfile);

data_voltagea = [(1:n-1)' angle(1 + X * conj(s))/pi*180 results.bus(PQnodes,VM) results.bus(PQnodes,VM)-margin_tightest/pi*180 results.bus(PQnodes,VM)+margin_tightest/pi*180];

save('-append', '-ascii', fname, 'data_voltagea');

