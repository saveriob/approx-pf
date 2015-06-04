clear all
%close all
clc

addpath('matpower4.1');
define_constants;

% GRID MODEL
mpc = loadcase('case_ieee123');

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

%e0 = zeros(n,1);
%e0(PCCindex) = 1;
%X2 = inv([L,e0;e0',0]); X2 = X2(1:n,1:n);
%X2 = X2(PQnodes,PQnodes);

%%%%%%%

s = mpc.bus(PQnodes,PD) + mpc.bus(PQnodes,GS) + 1j * (mpc.bus(PQnodes,QD) - mpc.bus(PQnodes,BS));

%%%%%%%%%

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = 1 + X * conj(s);
u_tild = abs(u_true - u_appr);


%%%%%%%%%

tb = [14 21];

range = -1.5:0.05:1.5;
[PP1 PP2] = meshgrid(range, range);
ZZ = zeros(size(PP1));

for p1 = 1:length(range)
for p2 = 1:length(range)

	mpc.bus(tb(1),PD) = range(p1);
	mpc.bus(tb(2),PD) = range(p2);

	results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));
	ZZ(p1,p2) = ...
		all(results.bus(PQnodes,VM) >= 0.95) ...
		&& all(results.bus(PQnodes,VM) <= 1.05);

end
end



plot(1:(n-1), results.bus(PQnodes,VM), 'k. ');
title('Voltage magnitude')
