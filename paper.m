%	This code generates the simulations included in
%
%	S. Bolognani, S. Zampieri (2014)
%	"On the existence and linear approximation of the power flow solution in power distribution networks"
%	Preprint available at http://arxiv.org/abs/1403.5031
%
%	This source code is distributed in the hope that it will be useful, but without any warranty.
%	We do request that publications in which this testbed is adopted, explicitly acknowledge that fact by citing the above mentioned paper.
%
%	MatLab OR GNU Octave, version 3.8.1 available at http://www.gnu.org/software/octave/
%	MATPOWER 4.1 available at http://www.pserc.cornell.edu/matpower/
%
%	tab width 4

clear all
close all
clc

more off

addpath('matpower5.1');
define_constants;

% Load case_ieee123, inspired by the IEEE 123 test feeder, with 
% - symmetric lines
% - balanced loads modeled as PQ buses
% - balanced shunt capacitors
% - switched in their normal position
% - ideal voltage regulators

% The modified testbed is distributed as case_ieee123 in the casefiles
% directory, and needs to be copied in the matpower directory.

mpc = loadcase('case_ieee123');

% Define useful constants

PCCindex = find(mpc.bus(:,BUS_TYPE)==3);
n = length(mpc.bus(:,BUS_TYPE));
PQnodes = setdiff(1:n,PCCindex);

% Build Laplacian L (neglecting shunt admittances)

nbr = size(mpc.branch,1);
nbu = size(mpc.bus,1);
L = zeros(nbu,nbu);

for br = 1:nbr
	br_F_BUS = mpc.branch(br,F_BUS); % FROM bus
	br_T_BUS = mpc.branch(br,T_BUS); % TO bus
	br_BR_R = mpc.branch(br,BR_R); % RESISTANCE
	br_BR_X = mpc.branch(br,BR_X); % INDUCTANCE
	br_Y = 1 / (br_BR_R + 1j * br_BR_X); % ADMITTANCE

	L(br_F_BUS, br_T_BUS) = br_Y;
	L(br_T_BUS, br_F_BUS) = br_Y;
	L(br_F_BUS, br_F_BUS) = L(br_F_BUS, br_F_BUS) - br_Y;
	L(br_T_BUS, br_T_BUS) = L(br_T_BUS, br_T_BUS) - br_Y;
end

% Build matrix X
X = inv(L(PQnodes,PQnodes));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figno = 1;

disp('FIGURE 1: Voltage magnitude');

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

s = mpc.bus(PQnodes,PD) + mpc.bus(PQnodes,GS) + 1j * (mpc.bus(PQnodes,QD) - mpc.bus(PQnodes,BS));

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);

figure(figno)
plot(1:(n-1), abs(u_true), 'ko ', 1:(n-1), abs(u_appr), 'k. ');
title('Voltage magnitude')
xlim([1 n-1]);


% check conditions of Theorem 1 for the existence of a practical solution

rho_2_star = max(arrayfun(@(idx) norm(X(idx,:),2), 1:size(X,1)));
fprintf(1,'X 2-star: %f\n', rho_2_star);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_2_star * norm(s,2));


% check if condition is verified with induced 2-norm

rho_2 = norm(X,2);
fprintf(1,'X 2-norm: %f\n', rho_2);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_2 * norm(s,2));


% check if condition is verified with p=1, q=Inf norms

rho_Inf_star = max(max(abs(X)));
fprintf(1,'X Inf-star: %f\n', rho_Inf_star);
fprintf(1,'s 1-norm: %f\n', norm(s,1));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_Inf_star * norm(s,1));


% plot error bands for different p-norms

ps = [1 2];
nps = length(ps);

margins = zeros(n-1,nps);

hold on

for p = 1:nps

	normq = arrayfun(@(idx) norm(X(idx,:),holder(ps(p))), 1:size(X,1));
	margins(:,p) = 4 * normq * max(normq) * norm(s,ps(p))^2;

	plot(1:(n-1), abs(u_appr)-margins(:,p), p, 1:(n-1), abs(u_appr)+margins(:,p), p);
	text(n, abs(u_appr(end))-margins(end,p), num2str(ps(p)));

end

ylim([0.5 1.1])

% Save voltage magnitude data

fname = 'data_paper_voltagem.data';

myfile=fopen(fname,'w');
fprintf(myfile,'bus appr real min1 max1 min2 max2');
fclose(myfile);

data_voltagem = [(1:n-1)'					...
		 abs(u_appr)						...
		 results.bus(PQnodes,VM)			...
		 abs(u_appr)-margins(:,1)			...
		 abs(u_appr)+margins(:,1)			...
		 abs(u_appr)-margins(:,2)			...
		 abs(u_appr)+margins(:,2)];

save('-append', '-ascii', fname, 'data_voltagem');

errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DC power flow comparison

figno = 2;

disp('FIGURE 2: Voltage angles');

figure(figno)

% plot 
% - true angles
% - approximate model
% - linear angles
% - DC power flow

plot(	1:(n-1), angle(u_true)/pi*180, 'ko ',		...
		1:(n-1), angle(u_appr)/pi*180, 'k. ',		...
		1:(n-1), imag(X * conj(s))/pi*180, 'kx ',	...
		1:(n-1), (imag(X) * real(s)) /pi*180, 'k+ ');
title('Voltage angle')
xlim([1 n-1]);


% Save voltage angle data

fname = 'data_paper_voltagea.data';

myfile=fopen(fname,'w');
fprintf(myfile,'bus appr real linear dc');
fclose(myfile);

data_voltagem = [(1:n-1)'					...
		 angle(u_appr)/pi*180				...
		 angle(u_true)/pi*180				...
		 imag(X * conj(s))/pi*180			...
		 (imag(X) * real(s)) /pi*180];

save('-append', '-ascii', fname, 'data_voltagem');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figno = 3;

disp('FIGURE 3: Voltage magnitude (uniform overload)');

mpc = loadcase('case_ieee123');

r = 2;

mpc.bus(PQnodes,PD) = r * mpc.bus(PQnodes,PD);
mpc.bus(PQnodes,GS) = r * mpc.bus(PQnodes,GS);
mpc.bus(PQnodes,QD) = r * mpc.bus(PQnodes,QD);
mpc.bus(PQnodes,BS) = r * mpc.bus(PQnodes,BS);

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

s = mpc.bus(PQnodes,PD) + mpc.bus(PQnodes,GS) + 1j * (mpc.bus(PQnodes,QD) - mpc.bus(PQnodes,BS));

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);

figure(figno)
plot(1:(n-1), abs(u_true), 'ko ', 1:(n-1), abs(u_appr), 'k. ');
title('Voltage magnitude')
xlim([1 n-1]);


% check conditions of Theorem 1 for the existence of a practical solution

rho_2_star = max(arrayfun(@(idx) norm(X(idx,:),2), 1:size(X,1)));

fprintf(1,'X 2-star: %f\n', rho_2_star);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_2_star * norm(s,2));


% check if condition is verified with induced 2-norm

rho_2 = norm(X,2);
fprintf(1,'X 2-norm: %f\n', rho_2);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_2 * norm(s,2));


% check if condition is verified with p=1, q=Inf norms

rho_Inf_star = max(max(abs(X)));
fprintf(1,'X Inf-star: %f\n', rho_Inf_star);
fprintf(1,'s 1-norm: %f\n', norm(s,1));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_Inf_star * norm(s,1));


% plot error bands for different p-norms

ps = [1 2];
nps = length(ps);

margins = zeros(n-1,nps);

hold on

for p = 1:nps

    normq = arrayfun(@(idx) norm(X(idx,:),holder(ps(p))), 1:size(X,1));
	margins(:,p) = 4 * normq * max(normq) * norm(s,ps(p))^2;

	plot(1:(n-1), abs(u_appr)-margins(:,p), p, 1:(n-1), abs(u_appr)+margins(:,p), p);
	text(n, abs(u_appr(end))-margins(end,p), num2str(ps(p)));

end

ylim([0.5 1.1])

% Save voltage magnitude data

fname = 'data_paper_voltagem_overload.data';

myfile=fopen(fname,'w');
fprintf(myfile,'bus appr real min1 max1 min2 max2');
fclose(myfile);

data_voltagem = [(1:n-1)'					...
		 abs(u_appr)						...
		 results.bus(PQnodes,VM)			...
		 abs(u_appr)-margins(:,1)			...
		 abs(u_appr)+margins(:,1)			...
		 abs(u_appr)-margins(:,2)			...
		 abs(u_appr)+margins(:,2)];

save('-append', '-ascii', fname, 'data_voltagem');

errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figno = 4;

disp('FIGURE 4: Voltage magnitude (spot overload)');

mpc = loadcase('case_ieee123');

r = 50;

mpc.bus(32,PD) = r * mpc.bus(32,PD);
mpc.bus(32,GS) = r * mpc.bus(32,GS);
mpc.bus(32,QD) = r * mpc.bus(32,QD);
mpc.bus(32,BS) = r * mpc.bus(32,BS);

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

s = mpc.bus(PQnodes,PD) + mpc.bus(PQnodes,GS) + 1j * (mpc.bus(PQnodes,QD) - mpc.bus(PQnodes,BS));

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);

figure(figno)
plot(1:(n-1), abs(u_true), 'ko ', 1:(n-1), abs(u_appr), 'k. ');
title('Voltage magnitude')
xlim([1 n-1]);


% check conditions of Theorem 1 for the existence of a practical solution

rho_2_star = max(arrayfun(@(idx) norm(X(idx,:),2), 1:size(X,1)));
fprintf(1,'X 2-star: %f\n', rho_2_star);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_2_star * norm(s,2));


% check if condition is verified with induced 2-norm

rho_2 = norm(X,2);
fprintf(1,'X 2-norm: %f\n', rho_2);
fprintf(1,'s 2-norm: %f\n', norm(s,2));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_2 * norm(s,2));


% check if condition is verified with p=1, q=Inf norms

rho_Inf_star = max(max(abs(X)));
fprintf(1,'X Inf-star: %f\n', rho_Inf_star);
fprintf(1,'s 1-norm: %f\n', norm(s,1));
fprintf(1,'4 * ||X|| * ||s|| = %f <? 1\n\n', rho_Inf_star * norm(s,1));


% plot error bands for different p-norms

ps = [1 2];
nps = length(ps);

margins = zeros(n-1,nps);

hold on

for p = 1:nps

    normq = arrayfun(@(idx) norm(X(idx,:),holder(ps(p))), 1:size(X,1));
	margins(:,p) = 4 * normq * max(normq) * norm(s,ps(p))^2;

	plot(1:(n-1), abs(u_appr)-margins(:,p), p, 1:(n-1), abs(u_appr)+margins(:,p), p);
	text(n, abs(u_appr(end))-margins(end,p), num2str(ps(p)));

end

ylim([0.5 1.1])


% Save voltage magnitude data

fname = 'data_paper_voltagem_spotoverload.data';

myfile=fopen(fname,'w');
fprintf(myfile,'bus appr real min1 max1 min2 max2');
fclose(myfile);

data_voltagem = [(1:n-1)'					...
		 abs(u_appr)						...
		 results.bus(PQnodes,VM)			...
		 abs(u_appr)-margins(:,1)			...
		 abs(u_appr)+margins(:,1)			...
		 abs(u_appr)-margins(:,2)			...
		 abs(u_appr)+margins(:,2)];

save('-append', '-ascii', fname, 'data_voltagem');

errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figno = 5;

disp('FIGURE 5: Voltage magnitude (PV bus)');

mpc = loadcase('case_ieee123');

PVbus = [15 51];

mpc.bus(PVbus,BUS_TYPE) = 2;
mpc.gen = 	[mpc.gen; ...
			[	PVbus(1)	0	0	200	-200	1	1	1	200	-200	0	0	0	0	0	0	0	0	0	0	0 ] ; ...
			[	PVbus(2)	0	0	200	-200	1	1	1	200	-200	0	0	0	0	0	0	0	0	0	0	0 ]];
mpc.gencost = [mpc.gencost; ...
				[ 2	0	0	3	0.01	40	0 ]; ...
				[ 2	0	0	3	0.01	40	0 ]];

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

s = mpc.bus(PQnodes,PD) + mpc.bus(PQnodes,GS) + 1j * (mpc.bus(PQnodes,QD) - mpc.bus(PQnodes,BS));

Q = -inv(imag(X(PVbus,PVbus)));
R = real(X(PVbus,PVbus)) * real(s(PVbus));
S = real(X(PVbus,setdiff(PQnodes,PVbus))) * real(s(setdiff(PQnodes,PVbus)));
T = imag(X(PVbus,setdiff(PQnodes,PVbus))) * imag(s(setdiff(PQnodes,PVbus)));

qPVbus = Q * (R + S + T);
s(PVbus) = s(PVbus) + 1j * qPVbus;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);

figure(figno)
plot(1:(n-1), abs(u_true), 'ko ', 1:(n-1), abs(u_appr), 'k. ');
title('Voltage magnitude')
xlim([1 n-1]);

ylim([0.9 1.05])

% Save voltage magnitude data

fname = 'data_paper_voltagem_pvbus.data';

myfile=fopen(fname,'w');
fprintf(myfile,'bus appr real');
fclose(myfile);

data_voltagem = [(1:n-1)'					...
		 abs(u_appr)						...
		 results.bus(PQnodes,VM)];

save('-append', '-ascii', fname, 'data_voltagem');

errorfigures;

