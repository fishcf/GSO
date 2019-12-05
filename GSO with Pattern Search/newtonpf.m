function [V, converged, i] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%NEWTONPF  Solves the power flow using a full Newton's method.
%   [V, converged, i] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. The bus voltage vector contains the set point for
%   generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. mpopt is a MATPOWER options vector which can be used to 
%   set the termination tolerance, maximum number of iterations, and 
%   output options (see 'help mpoption' for details). Uses default
%   options if this parameter is not given. Returns the final complex
%   voltages, a flag which indicates whether it converged or not, and
%   the number of iterations performed.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/24/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% default arguments
if nargin < 7
	mpopt = mpoption;
end

%% options
tol		= mpopt(2);
max_it	= mpopt(3);
verbose	= mpopt(31);

%% initialize
j = sqrt(-1);
converged = 0;
i = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

%% set up indexing for updating V
npv	= length(pv);
npq	= length(pq);
j1 = 1;			j2 = npv;			%% j1:j2 - V angle of pv buses
j3 = j2 + 1;	j4 = j2 + npq;		%% j3:j4 - V angle of pq buses
j5 = j4 + 1;	j6 = j4 + npq;		%% j5:j6 - V mag of pq buses

%% evaluate F(x0)
mis = V .* conj(Ybus * V) - Sbus;
F = [	real(mis([pv; pq]));
		imag(mis(pq))	];

%% check tolerance
normF = norm(F, inf);
if verbose > 1
	fprintf('\n it    max P & Q mismatch (p.u.)');
	fprintf('\n----  ---------------------------');
	fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
	converged = 1;
	if verbose > 1
		fprintf('\nConverged!\n');
	end
end

%% do Newton iterations
while (~converged & i < max_it)
	%% update iteration counter
	i = i + 1;
	
	%% evaluate Jacobian
	[dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);
	
	%% selecting a subset of rows of a large sparse matrix is very slow
	%% in Matlab 5 (but not Matlab 4 ... go figure), but selecting a
	%% subset of the columns is fast, and so is transposing, so instead
	%% of doing this ...
% 	j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
% 	j12 = real(dSbus_dVm([pv; pq], pq));
% 	j21 = imag(dSbus_dVa(pq, [pv; pq]));
% 	j22 = imag(dSbus_dVm(pq, pq));

	%% ... we do the equivalent thing using
	%% a temporary matrix and transposing
	temp = real(dSbus_dVa(:, [pv; pq]))';
	j11 = temp(:, [pv; pq])';
	temp = real(dSbus_dVm(:, pq))';
	j12 = temp(:, [pv; pq])';
	temp = imag(dSbus_dVa(:, [pv; pq]))';
	j21 = temp(:, pq)';
	temp = imag(dSbus_dVm(:, pq))';
	j22 = temp(:, pq)';
	
	J = [	j11 j12;
			j21 j22;	];

	%% compute update step
	dx = -(J \ F);

	%% update voltage
	Va(pv) = Va(pv) + dx(j1:j2);
	Va(pq) = Va(pq) + dx(j3:j4);
	Vm(pq) = Vm(pq) + dx(j5:j6);
	V = Vm .* exp(j * Va);

	%% evalute F(x)
	mis = V .* conj(Ybus * V) - Sbus;
	F = [	real(mis(pv));
			real(mis(pq));
			imag(mis(pq))	];

	%% check for convergence
	normF = norm(F, inf);
	if verbose > 1
	fprintf('\n%3d        %10.3e', i, normF);
	end
	if normF < tol
		converged = 1;
		if verbose
			%fprintf('\nNewton''s method power flow converged in %d iterations.\n', i);
		end
	end
end

if verbose
	if ~converged
		fprintf('\nNewton''s method power did not converge in %d iterations.\n', i);
	end
end
