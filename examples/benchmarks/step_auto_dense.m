% Benchmark for dense matrix propagation backends. Compares the
% default Taylor series with the heuristic selection of faster
% exponential-action algorithms. Syntax:
%
%                  step_auto_dense()
%
% a.acharya@soton.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=step_auto_dense.m>

function step_auto_dense()

% Set problem dimensions
dims=[128,256,512];

% Set matrix norms
norms=[2,10,60];

% Bootstrap the spin system
spin_system=bootstrap('hush');

% Loop over dimensions
for n=1:numel(dims)

    % Loop over norms
    for k=1:numel(norms)

        % Generate a random dense Hamiltonian
        H=randn(dims(n),dims(n))+1i*randn(dims(n),dims(n));

        % Make the Hamiltonian Hermitian
        H=(H+H')/2;

        % Scale the Hamiltonian
        H=norms(k)*H/norm(H,1);

        % Generate a random initial state
        rho=randn(dims(n),1)+1i*randn(dims(n),1);

        % Normalise the initial state
        rho=rho/norm(rho,2);

        % Set time step
        dt=1.0;

        % Set backend to default
        spin_system.sys.expmv_backend='default';

        % Time the default backend
        tic;rho_def=step(spin_system,H,rho,dt);t_def=toc;

        % Set backend to auto
        spin_system.sys.expmv_backend='auto';

        % Time the auto backend
        tic;rho_auto=step(spin_system,H,rho,dt);t_auto=toc;

        % Compute the error
        err=norm(rho_def-rho_auto,2)/norm(rho_def,2);

        % Report the results
        fprintf('Dim: %4d, Norm: %2d, Default: %6.4fs, Auto: %6.4fs, Speedup: %5.2fx, Error: %1.2e\n',dims(n),norms(k),t_def,t_auto,t_def/t_auto,err);

    end

end

end

% "The most important thing in the world is to be able to
%  be what you are."
%
%  Harlan Ellison
