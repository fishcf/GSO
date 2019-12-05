function [Fuelcost,valve_value,Emission,TaxCO2,Totalloss,Lmax, TotalPg,SumVoltageDeviations,Vlim,Qglim,Slim,BusV,gen,isConverge] = calpfnew(casename,Individual)
global bus;
% [Fuelcost,SumVoltageDeviations,Lmax,Vlim,Qglim] = calpfnew(casename,Individual)
%                          
%   Run power flow with a set of control varibles. 
%   
%Input Arguments:
%   casename       - the name of the evaluation .m function
%   Individual      - a set of control varibles from algorithms

% Modified power flow program for Matlab  
% Copyright (C) 2003 Shan He, the University of Liverpool
% Intelligence Engineering & Automation Group

%RUNPF  Runs a power flow.
%
%   [baseMVA, bus, gen, branch, success, et] = runpf(casename, mpopt, fname)
%
%   Runs a full Newton's method power flow where casename is the name of
%   the m-file (without the .m extension) containing the power flow data,
%   and mpopt is a MATPOWER options vector (see 'help mpoption' for details).
%   Uses default options if 2nd parameter is not given, and 'case' if 1st
%   parameter is not given. The results may optionally be printed to a file
%   (appended if the file exists) whose name is given in fname (in addition
%   to printing to STDOUT). Optionally returns the final values of baseMVA,
%   bus, gen, branch, success, and et.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/24/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.
tic;
%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
%% default arguments
%% define the shunt VAR compensation buses indices
%  QcIndices = [10 12 15 17 20 21 23 24 29];


fname = '';
mpopt = mpoption;		%% use default options


%% options
alg = 1;

%% read data & convert to internal bus numbering
[baseMVA, bus, gen, branch, area, gencost] = feval(casename);

xfmr = find(branch(:, TAP));
genbus = [find(bus(:,BUS_TYPE)==3); find(bus(:,BUS_TYPE)==2)];
QcIndices = find(bus(:,BS));

% Decode parameter from individual
cutpoint=size(genbus,1)-1;
Pg = Individual(1:cutpoint);
cutpoint2 = cutpoint+size(bus(genbus,VMAX),1);
cutpoint = cutpoint+1;
Vg = Individual(cutpoint:cutpoint2);
cutpoint = cutpoint2+size(xfmr,1);
cutpoint2=cutpoint2+1;
% T = .9+(floor(Individual(cutpoint2:cutpoint))-1)*.025;
% %T = Individual(cutpoint2:cutpoint);


T = .9+(floor(Individual(cutpoint2:cutpoint))-1)*.0125;  %……17个离散值



cutpoint2=cutpoint+size(QcIndices,1);
% For some cases, there is no shunt at all.
if cutpoint~=cutpoint2
    cutpoint = cutpoint+1;
%     Qc = floor(Individual(cutpoint:cutpoint2));
    Qc = (floor(Individual(cutpoint:cutpoint2)));
    
    if(Qc==6)
        Qc=5;
    end
    bus(QcIndices,BS)=Qc;
end


gen(2:size(genbus,1),PG) = Pg;
gen(:,VG) = Vg;
branch(find(branch(:,TAP)),TAP) = T;

[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% build admittance matrices
%   Yf and Yt, when multiplied by a complex voltage vector, yield the vector
%   currents injected into each line from the "from" and "to" buses
%   respectively of each line.
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% compute complex bus power injections (generation - load)
%returns the vector of complex bus power injections
Sbus = makeSbus(baseMVA, bus, gen);

%% generator info
on = find(gen(:, GEN_STATUS));				%% which generators are on?


%% initialize V and Pg from data from case file
%发电机为PV节点
% V0	= ones(size(bus, 1), 1);		%% flat start
V0	= bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));%计算节点复电压,但为何一同计算PQ节点？
%更新发电机所在节点的电压,y=V*V0/|V0|,都是标幺值,为何不一样
V0(gen(on, GEN_BUS)) = gen(on, VG) ./ abs(V0(gen(on, GEN_BUS))).* V0(gen(on, GEN_BUS));


%% run the power flow
t0 = clock;
if alg == 1
	[V, success, iterations] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt);
elseif alg == 2 | alg == 3
	[Bp, Bpp] = makeB(baseMVA, bus, branch, alg);
	[V, success, iterations] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt);
else
	error('Only Newton''s method and fast-decoupled power flow algorithms currently implemented.');
end
%gen(:,PG)
%% compute flows etc.
[bus, gen, branch] = pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq);
et = etime(clock, t0);

%%% convert back to original bus numbering & print results
[bus, gen, branch] = int2ext(i2e, bus, gen, branch);

LoadBusV=bus(find(bus(:,BUS_TYPE)==1),VM);
BusV=bus(:,VM);
SumVoltageDeviations=sum(abs(LoadBusV-1));

Vmin=bus([find(bus(:,BUS_TYPE)==1);find(bus(:,BUS_TYPE)==3)],VMIN);
Vmax=bus([find(bus(:,BUS_TYPE)==1);find(bus(:,BUS_TYPE)==3)],VMAX);
Qmax=gen(:,QMAX);
Qmin=gen(:,QMIN);

Pmax=gen(:,PMAX);
Pmin=gen(:,PMIN);
%Qg=gen(:,QG)

%平衡节点和PQ节点
ConstrainedBusV=[bus(find(bus(:,BUS_TYPE)==3),VM); bus(find(bus(:,BUS_TYPE)==1),VM)];
%衡量超过最值约束的电压值
Vlim=(abs(ConstrainedBusV)>Vmax).*(abs(ConstrainedBusV)-Vmax).^2 + (abs(ConstrainedBusV)<Vmin).*(abs(ConstrainedBusV)-Vmin).^2;
%abs(LoadBusV);
Qglim=(gen(:,QG)>Qmax).*(gen(:,QG)-Qmax).^2 + (gen(:,QG)<Qmin).*(gen(:,QG)-Qmin).^2;

nb = size(bus, 1); %% number of buses
ng = size(gen, 1);  %% number of generators

V = bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
Si=zeros(nb,1);  
Yjj=zeros(nb,1);
Sicorr=zeros(nb,1);
Zbus=Ybus^-1;



% find out loadbus indices
%找出所有PQ节点
a=0;
for i = 1:nb
    if bus(i, BUS_TYPE)==1
        a=a+1;
	    LoadIndic(a,1)=i;
        %Ynewbus(a,a)=Ybus(i,i);        
    end
end


% construct YLL submatrix
%找出所有PQ节点之间的导纳,24*24
YLLbus=zeros(a,a);
for i=1:a
    for j=1:a
        YLLbus(i,j)=Ybus(LoadIndic(i,1), LoadIndic(j,1));
    end
end


% find out generator indices
b=0;
for i = 1:nb
    if bus(i, BUS_TYPE)==2 | bus(i, BUS_TYPE)==3
        b=b+1;
        genbusnum(b,1)=i;        
    end
end

 % construct YLG submatrix
 %寻找PQ节点与发电机节点间的导纳
YLGbus=zeros(a,ng);
for i=1:a
    for j=1:b
        YLGbus(i,j)=Ybus(LoadIndic(i,1),genbusnum(j,1));
    end
end

% construct submatrix F of hybird matrix H
Fbus=-YLLbus^-1*YLGbus;

% calculation of L-index
FtimesV=zeros(a,1);
for j=1:a   % j is the load bus
    for i=1:ng  % i is the generator bus
        FtimesV(j,1)=FtimesV(j,1)+Fbus(j,i)*V(genbusnum(i,1),1);
    end
    if bus(LoadIndic(j,1),PD)~=0 & bus(LoadIndic(j,1),PQ)~=0 
        Lindex(j,1)=abs(1-FtimesV(j,1)/V(LoadIndic(j,1),1)); 
    end

end

Lindex;
Lmax=max(Lindex);%L是？？？



% 
% e2i = zeros(max(bus(:, BUS_I)), 1);		%% need internal bus numbering for a second
% e2i(bus(:, BUS_I)) = [1:size(bus, 1)]';
% 
% nl = size(branch, 1);	%% number of branches
% tap = ones(nl, 1);								%% default tap ratio = 1 for lines
% xfmr = find(branch(:, TAP));					%% indices of transformers
% tap(xfmr) = branch(xfmr, TAP);					%% include transformer tap ratios
% tap = tap .* exp(j*pi/180 * branch(:, SHIFT));	%% add phase shifters
% loss = baseMVA * (abs(V(e2i(branch(:, F_BUS))) ./ tap - V(e2i(branch(:, T_BUS)))) .^ 2 ./ ...
% 				(branch(:, BR_R) - j * branch(:, BR_X)));
% loss =abs(loss);
% totalloss=sum(real(loss));
TotalPg=sum(gen(:,PG));
TotalPd=sum(bus(:,3));
Totalloss=TotalPg-TotalPd;

% gen(:,PG)=gen(:,PG)./baseMVA;


%     [Pg,Vg,T,Qc] = DecodeParameters(x);
%     gen(:,PG);
%     SumVoltageDeviations;
%     TotalFuel = sum( b'.*gen(:,PG)/baseMVA + c'.*(gen(:,PG)/baseMVA).^2) +0*SumVoltageDeviations
%     gencost;
%     
    
    Fuelcost = sum(totcost(gencost, gen(:,PG)));
     gen(:,PG);
     emission=[1 22.983 -1.1000 0.0126
            2 25.313 -0.1000 0.0200
            3 25.505 -0.0100 0.0270
            4 24.900 -0.0050 0.0291
            5 24.700 -0.0040 0.0290
            6 25.300 -0.0055 0.0271];
      valve= [18 0.037 
              16 0.038 
              14 0.04 
              12 0.045 
              13 0.042 
              13.5 0.041];
       Pmin= [50.0000; 20.0000; 15.0000; 10.0000; 10.0000; 12.0000];
       CO2=[40  0.2  0.00004 
            50  0.3  0.00005 
            80  0.12  0.000024 
            2462.4  48  0.0084 
            2500  50  0.009 
            1.248  0.234  0.00003432];
        CO2factor=[3.1604  3.1604  3.1604  1.84 1.84 2.8523]';
        TaxCO2=0.0095*sum(CO2factor.*(CO2(:,1)+CO2(:,2).*gen(:,PG)+ CO2(:,3).*(gen(:,PG).^2)));
%         emission(:,3).*gen(:,PG)
%         sum(emission(:,2))
%         sum(emission(:,3).*gen(:,PG))
%         sum(emission(:,4).*(gen(:,PG).^2))
      valve_value=sum(abs(valve(:,1).*sin(valve(:,2).*(Pmin-gen(:,PG)))));
      Emission= sum(emission(:,2))+ sum(emission(:,3).*gen(:,PG))+sum(emission(:,4).*(gen(:,PG).^2));
%% compute branch power flows
	Sf = V(branch(:, F_BUS)) .* conj(Yf * V);	%% complex power injected at "from" bus (p.u.)
	St = V(branch(:, T_BUS)) .* conj(Yt * V);	%% complex power injected at "to" bus (p.u.)

%      a=branch(:, RATE_A) / baseMVA;
%      b=abs(Sf);
%% inequality constraints (line flow limits)
% 	Sflim = (abs(Sf) - branch(:, RATE_A) / baseMVA);	%% branch apparent power limits (from bus)
% 	Stlim =	(abs(St) - branch(:, RATE_A) / baseMVA);	%% branch apparent power limits (to bus)
    Slim=(abs(Sf-St)-branch(:, RATE_A) / baseMVA);
    
if fname
	[fd, msg] = fopen(fname, 'at');
	if fd == -1
		error(msg);
	else
		printpf(baseMVA, bus, gen, branch, [], success, et, fd, mpopt);
		fclose(fd);
	end
end

%% 判断是否收敛
isConverge=1;
isSlim=max(Slim)>0.2;
% isVlim=(any(abs(ConstrainedBusV)>Vmax))||(any(abs(ConstrainedBusV)<Vmin));
% isQglim=any(gen(:,QG)>Qmax) || any(gen(:,QG)<Qmin);
isPglim=any(gen(:,PG)>Pmax) || any(gen(:,PG)<Pmin);
% if(~success||isSlim||isVlim||isQglim||isPglim)
if(~success||isSlim||isPglim)
    isConverge=0;
end

%printpf(baseMVA, bus, gen, branch, [], success, et, 1, mpopt);
%gen(:,PG)
%% this is just to prevent it from printing baseMVA
%% when called with no output arguments
%bus
if nargout, MVAbase = baseMVA; end

return;
