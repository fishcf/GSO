function [cur,Slim_max]=CalFitnessnew(casename,x,numofscens)
global bus;
global  windenvironment; 
% [Fitness,SumVoltageDeviations,Vlim,Qglim]=CalFitnessnew(casename,x)
%                          
%   Run power flow with a set of control varibles. 
%   
% Input Arguments:
%   casename       - the name of the evaluation .m function
%   x      - a set of control varibles from algorithms

% Modified power flow program for Matlab  
% Copyright (C) 2003 Shan He, the University of Liverpool
% Intelligence Engineering & Automation Group


%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, ~, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;


lambda1=100;
lambda2=0.01;
global beta;


%x = InitializePop(Bound, 1)
%Bound=[20 80; 15 50; 10 35; 10 30; 12 40; [.95*ones(6,1) 1.1*ones(6,1)]; [1*ones(4,1) 9*ones(4,1)]; [0*ones(9,1) 5*ones(9,1)]];
%b=[200 175 100 325 300 300];
%c=[37.5 175 625 83.4 250 250];

if nargin==1
    [baseMVA, bus, gen, branch, area, gencost] = feval(casename);
    xfmr = find(branch(:, TAP));
    genbus = [find(bus(:,BUS_TYPE)==3); find(bus(:,BUS_TYPE)==2)];
    QcIndices = find(bus(:,BS));

    %  Vmax and Vmin
    Bound=[gen(2:size(gen,1),PMIN) gen(2:size(gen,1),PMAX); bus(genbus,VMIN) bus(genbus,VMAX); ... 
      % Tap position                                    % Shunt
        1*ones(size(xfmr,1),1) 17*ones(size(xfmr,1),1); 0*ones(size(QcIndices,1),1) 6*ones(size(QcIndices,1),1)];
    
else
    [Fuelcost,valve_value,Emission,TaxCO2,Totalloss,Lmax, TotalPg,SumVoltageDeviations,Vlim,Qglim,Slim,BusV,gen,isConverge] = calpfnew(casename,x);
    
%      
     %计算最大超载率
     Slim_max=max(Slim);
     
     windpositions=[2 7 10 16 24];
     cur=0;

     counts=0;
     lambda1=100;
     lambda2=0.01;
     lambda1=lambda1*510;
     lambda2=lambda2*100;
     cur= 1*Fuelcost + 1*lambda1*sum(Vlim)+1*lambda2*sum(Qglim);
     %apply curtailment
     while(~isConverge)   
          counts=counts+1;
          bus(windpositions,3) = bus(windpositions,3)+windenvironment(numofscens,:)'*0.025;
% 
%          cur=cur+sum(windenvironment(numofscens,:)')*0.025;
         [Fuelcost,valve_value,Emission,TaxCO2,Totalloss,Lmax, TotalPg,SumVoltageDeviations,Vlim,Qglim,Slim,BusV,gen,isConverge] = calpfnew(casename,x);
    
         Slim_max=max(Slim);
         if(any(counts>=40))
             break;
         end
     end
     
end

%if cur
