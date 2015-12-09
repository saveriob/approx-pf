%	This code generates a extra set of simulations of the model presented in
%
%	S. Bolognani, S. Zampieri
%	"On the existence and linear approximation of the power flow solution in power distribution networks"
%	to appear on IEEE Transactions on Power Systems.
%	doi: 10.1109/TPWRS.2015.2395452
%	Preprint available at http://arxiv.org/abs/1403.5031
%
%	This source code is distributed in the hope that it will be useful, but without any warranty.
%	We do request that publications in which this testbed is adopted, explicitly acknowledge that fact by citing the above mentioned paper.
%
%	MatLab OR GNU Octave, version 3.8.1 available at http://www.gnu.org/software/octave/
%	MATPOWER 5.1 available at http://www.pserc.cornell.edu/matpower/
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

% As described in the attached technical note, the following simulations will be performed:
% 1) nominal scenario
% 2) uniform overload
% 3) lumped overload
% 4) shunt capacitors
% 5) voltage regulation
% 6) tap changer
% 7) PV buses
% 8) high R/X ratio

% For all these cases, voltage magnitude and voltage phases will be plotted.

% Define useful constants

PCCindex = find(mpc.bus(:,BUS_TYPE)==3);
n = length(mpc.bus(:,BUS_TYPE));
PQnodes = setdiff(1:n,PCCindex);

% Build Laplacian L (neglecting shunt admittances)

nbr = size(mpc.branch,1);
nbu = size(mpc.bus,1);
L = zeros(nbu,nbu);

for br = 1:nbr
	br_F_BUS = mpc.branch(br,F_BUS);
	br_T_BUS = mpc.branch(br,T_BUS);
	br_BR_R = mpc.branch(br,BR_R);
	br_BR_X = mpc.branch(br,BR_X);
	br_Y = 1 / (br_BR_R + 1j * br_BR_X);

	L(br_F_BUS, br_T_BUS) = -br_Y;
	L(br_T_BUS, br_F_BUS) = -br_Y;
	L(br_F_BUS, br_F_BUS) = L(br_F_BUS, br_F_BUS) + br_Y;
	L(br_T_BUS, br_T_BUS) = L(br_T_BUS, br_T_BUS) + br_Y;
end

% Build matrix X

X = inv(L(PQnodes,PQnodes));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 1

figno = 1;

disp('FIGURE 1: Nominal scenario');

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = 1;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = UPCC + X * conj(s);

u_true_nom = u_true;
u_appr_nom = u_appr;

figure(figno)
subplot(211)
	plot(1:(n-1), abs(u_true), 'ko ', 1:(n-1), abs(u_appr), 'k. ');
	title('Voltage magnitude')
	xlim([1 n-1]);
subplot(212)
	plot(1:(n-1), angle(u_true), 'ko ', 1:(n-1), angle(u_appr), 'k. ');
	title('Voltage angle')
	xlim([1 n-1]);

savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 2

figno = 2;

disp('FIGURE 2: Uniform overload');

mpc = loadcase('case_ieee123');

r = 2;

mpc.bus(PQnodes,PD) = r * mpc.bus(PQnodes,PD);
mpc.bus(PQnodes,GS) = r * mpc.bus(PQnodes,GS);
mpc.bus(PQnodes,QD) = r * mpc.bus(PQnodes,QD);
mpc.bus(PQnodes,BS) = r * mpc.bus(PQnodes,BS);

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = 1;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = UPCC + X * conj(s);

plotdata;
savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 3

figno = 3;

disp('FIGURE 3: Lumped overload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

mpc = loadcase('case_ieee123');

r = 50;

mpc.bus(32,PD) = r * mpc.bus(32,PD);
mpc.bus(32,GS) = r * mpc.bus(32,GS);
mpc.bus(32,QD) = r * mpc.bus(32,QD);
mpc.bus(32,BS) = r * mpc.bus(32,BS);

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = 1;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = UPCC + X * conj(s);

plotdata;
savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 4

figno = 4;

disp('FIGURE 4: Shunt capacitor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

mpc = loadcase('case_ieee123');

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);

mpc.bus(26,BS) = 0.600;
mpc.bus(28,BS) = 0.050;
mpc.bus(29,BS) = 0.050;
mpc.bus(30,BS) = 0.050;

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

UPCC = 1;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);

L4 = L;
L4(26,26) = L4(26,26) + 1j * 0.600;
L4(28,28) = L4(28,28) + 1j * 0.050;
L4(29,29) = L4(29,29) + 1j * 0.050;
L4(30,30) = L4(30,30) + 1j * 0.050;

X4 = inv(L4(PQnodes,PQnodes));
w = -X4 * L4(PQnodes,PCCindex);

u_appr = UPCC*w + X4 * inv(diag(conj(w))) * conj(s);

plotdata;
savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 5

figno = 5;

disp('FIGURE 5: Voltage regulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

mpc = loadcase('case_ieee123');

t = 120/124;
mpc.branch(17,TAP) = t;

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = 1;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);

% Build grid matrices for tap position t

reg_buses = [];
for bus = 1:n-1
	if (abs(X(bus,17)-X(17,17))<1e-5)
		reg_buses = [reg_buses bus];
	end
end

L5 = zeros(nbu,nbu);

for br = 1:nbr
	br_F_BUS = mpc.branch(br,F_BUS);
	br_T_BUS = mpc.branch(br,T_BUS);
	br_BR_R = mpc.branch(br,BR_R);
	br_BR_X = mpc.branch(br,BR_X);
	br_Y = 1 / (br_BR_R + 1j * br_BR_X);

	if (ismember(br_F_BUS, reg_buses) || ismember(br_T_BUS, reg_buses))
		tt = t^2;
	else
		tt = 1;
	end

	L5(br_F_BUS, br_T_BUS) = - tt * br_Y;
	L5(br_T_BUS, br_F_BUS) = - tt * br_Y;
	L5(br_F_BUS, br_F_BUS) = L5(br_F_BUS, br_F_BUS) + tt * br_Y;
	L5(br_T_BUS, br_T_BUS) = L5(br_T_BUS, br_T_BUS) + tt * br_Y;
end

X5 = inv(L5(PQnodes,PQnodes));

d = ones(n-1,1);
d(reg_buses) = 1/t;

u_appr =  diag(d) * (UPCC + X5 * conj(s));

plotdata;
savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 6

figno = 6;

disp('FIGURE 6: Tap changer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

mpc = loadcase('case_ieee123');

tp = 124/120;

mpc.gen(1,VG) = tp;

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = tp;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = UPCC + X/tp * conj(s);

plotdata;
savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 7

figno = 7;

disp('FIGURE 7: PV buses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

mpc = loadcase('case_ieee123');

PVbus = [15 51];

mpc.bus(PVbus,BUS_TYPE) = 2;
mpc.gen = 	[mpc.gen; ...
			PVbus'	ones(size(PVbus'))*[0	0	200	-200	1	1	1	200	-200	0	0	0	0	0	0	0	0	0	0	0 ]];
mpc.gencost = 	[mpc.gencost; ...
				ones(size(PVbus'))*[ 2	0	0	3	0.01	40	0 ]];

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = 1;

Q = -inv(imag(X(PVbus,PVbus)));
R = real(X(PVbus,PVbus)) * real(s(PVbus));
S = real(X(PVbus,setdiff(PQnodes,PVbus))) * real(s(setdiff(PQnodes,PVbus)));
T = imag(X(PVbus,setdiff(PQnodes,PVbus))) * imag(s(setdiff(PQnodes,PVbus)));

qPVbus = Q * (R + S + T);
s(PVbus) = s(PVbus) + 1j * qPVbus;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = UPCC + X * conj(s);

plotdata;
savedata;
errorfigures;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CASE 8

figno = 8;

disp('FIGURE 8: High R/X ratio %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

mpc = loadcase('case_ieee123_ug');

L8 = zeros(nbu,nbu);

for br = 1:nbr
	br_F_BUS = mpc.branch(br,F_BUS);
	br_T_BUS = mpc.branch(br,T_BUS);
	br_BR_R = mpc.branch(br,BR_R);
	br_BR_X = mpc.branch(br,BR_X);
	br_Y = 1 / (br_BR_R + 1j * br_BR_X);

	L8(br_F_BUS, br_T_BUS) = -br_Y;
	L8(br_T_BUS, br_F_BUS) = -br_Y;
	L8(br_F_BUS, br_F_BUS) = L8(br_F_BUS, br_F_BUS) + br_Y;
	L8(br_T_BUS, br_T_BUS) = L8(br_T_BUS, br_T_BUS) + br_Y;
end

% Build matrix X

X8 = inv(L8(PQnodes,PQnodes));

results = runpf(mpc, mpoption('VERBOSE', 0, 'OUT_ALL',0));

full_s = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen);
s = full_s(PQnodes);
UPCC = 1;

u_true = results.bus(PQnodes,VM) .* exp(1j * results.bus(PQnodes,VA)/180*pi);
u_appr = UPCC + X8 * conj(s);

plotdata;
savedata;
errorfigures;

