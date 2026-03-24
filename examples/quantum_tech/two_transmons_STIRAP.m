% Basic single transmon system with Duffing model inter-
% actions, parameters and model based off of (https://doi.org/10.1038/ncomms10628)
% by K.S.Kumar et al. - STIRAP PROTOCOL
%
% c.musselwhite@soton.ac.uk
%
% STILL VERY MUCH A WIP!
% 03/03/2026 ADJUSTMENT - Shift frequencies in terms of MHz (1e6) and times in
% terms of μs (1e-6).

function two_transmons_STIRAP()

% Magnet field
sys.magnet=0;

% Particle specification
sys.isotopes={'T3'};

% Formalism and basis
bas.formalism='zeeman-liouv';
bas.approximation='none';

% Spinach housekeeping
spin_system=create(sys);
spin_system=basis(spin_system,bas);
%spin_system=assume(spin_system,'duffing');

% Get elementary operators
%CrA=operator(spin_system,'C',1);
%AnA=operator(spin_system,'A',1);

% Hamiltonian parameters
%deltas=[10e6, -10e6];                       % Two-photon resonance is assumed - LOOK INTO THIS!
deltas=[10, -10];
%omegas=2*pi*[5.27e9 4.82e9];
%alphas=2*pi*-450e6;

% Builds Gaussian Pulse
% ORIGINAL VALUES
%t=1e-9*linspace(-150,150,2000);
%ts=-90e-9;
%sigma=45e-9;
%omega01=2*pi*43.4e6;
%omega12=2*pi*38.2e6;
% SCALED VALUES
t=1e-3*linspace(-150,150,300);
ts=-90e-3;
sigma=45e-3;
omega01=2*pi*43.4;
omega12=2*pi*38.2;
omega01_t=omega01*exp(-((t.^2)/(2*sigma^2)));
omega12_t=omega12*exp(-(((t-(ts/2)).^2)/(2*sigma^2)));

% Build the Drift Hamiltonian -> ASK ILYA ABOUT CONTROL HAMILTONIAN
%H=    deltas*operator(spin_system,'N',1)+...
%  (alphas/2)*operator(spin_system,'CCAA',1);

H=sparse([0 0 0; 0 deltas(1) 0; 0 0 (deltas(1)+deltas(2))]);    % NO DETUNING CASE - ALL ELEMENTS OF DRIFT = 0 (definitely incorrect)

% Build control operators
H01=sparse([0 1 0; 1 0 0; 0 0 0]/2);
H12=sparse([0 0 0; 0 0 1; 0 1 0]/2);

%C_x=(CrA+AnA)/2; 
%C_y=(CrA-AnA)/(2i);

% Build offset operators
%O_A=operator(spin_system,'N',1);     % Placeholder for detuning offset operator

% Build source and destination states
rho_init = state(spin_system,'BL1',1);   % |0>
rho_targ = state(spin_system,'BL3',1);   % |2>

rho_init=rho_init/norm(rho_init,'fro');
rho_targ=rho_targ/norm(rho_targ,'fro');

% Unit fidelity is Sorensen bound
%rho_targ=rho_targ/sorensen(rho_init,rho_targ);

% Define control parameters
%control.drifts={{H}};
control.drifts={{hilb2liouv(H,'comm')}};                             % Drift
control.operators={hilb2liouv(H01,'comm'),hilb2liouv(H12,'comm')};
%control.operators={C_x,C_y};                      % Controls
%control.off_ops={O_A};                        % Offset operator
%control.offsets={linspace(-0.5e6,0.5e6,5)};
control.rho_init={rho_init};                      % Starting state
control.rho_targ={rho_targ};                      % Destination state
control.pwr_levels=2*pi*[30 35 40 45 50];       % Pulse power ensemble
%control.pwr_levels=2*pi*[40 45 50 55 60]*1e6;
control.pulse_dt=1e-5*ones(1,300);                % Slice durations - SCALED
%control.pulse_dt=1e-13*ones(1,2000);               % Slice durations - ORIGINAL
control.penalties={'NS'};                         % Penalties
control.p_weights=10e-7;                            % Penalty weights
control.method='goodwin';                           % Optimisation method
control.max_iter=100;                             % Termination condition
control.parallel='ensemble';                      % Parallelisation mode

% Plots during optimisation
control.plotting={'xy_controls','spectrogram','robustness'};

% Spinach housekeeping
spin_system=optimcon(spin_system,control);

% Initial guess - random
pulse=[omega01_t; omega12_t];

% Run the optimisation, get normalised pulse
fmaxnewton(spin_system,@grape_xy,pulse);

end