function [ref, pv, pq] = bustypes(bus, gen)
%BUSTYPES   Builds lists of each type of bus (ref, pv, pq).
%   [ref, pv, pq] = bustypes(bus, gen)
%   Generators with "out-of-service" status are treated as PQ buses with
%   zero generation (regardless of Pg/Qg values in gen).

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    9/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% get generator status
bus_gen_status = zeros(size(bus, 1), 1);
bus_gen_status(gen(:, GEN_BUS)) = gen(:, GEN_STATUS);

%% form index lists for slack, PV, and PQ buses
ref	= find(bus(:, BUS_TYPE) == REF & bus_gen_status);	%% reference bus index
pv	= find(bus(:, BUS_TYPE) == PV  & bus_gen_status);	%% PV bus indices
pq	= find(bus(:, BUS_TYPE) == PQ | ~bus_gen_status);	%% PQ bus indices

%% pick a new reference bus if for some reason there is none (may have been shut down)
if isempty(ref)
	ref = pv(1);				%% use the first PV bus
	pv = pv(2:length(pv));		%% take it off PV list
end

return;
