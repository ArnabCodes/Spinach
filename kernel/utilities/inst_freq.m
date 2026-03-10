% Instantaneous frequency trajectory from a complex time-domain signal.
% Uses product-phase finite differences, syntax:
%
%                    freq=inst_freq(signal,dt)
%
% Parameters:
%
%    signal - complex row or column vector with time-domain signal
%
%    dt     - time step duration between signal points (seconds)
%
% Outputs:
%
%    freq   - instantaneous frequency trajectory (Hz), same size as signal
%
% Note: central second-order phase-difference approximation is used in the
%       interior, and second-order one-sided approximations are used at the
%       boundaries; all formulas use one-step phase differences
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=inst_freq.m>

function freq=inst_freq(signal,dt)

% Check consistency
grumble(signal,dt);

% Convert input into a column vector
signal_col=signal(:);

% Preallocate output array
freq=zeros(size(signal_col));

% Compute interior points using central second-order formula
% from one-step phase differences
phase_fwd=angle(signal_col(3:end).*conj(signal_col(2:end-1)));
phase_bwd=angle(signal_col(2:end-1).*conj(signal_col(1:end-2)));
freq(2:end-1)=(phase_fwd+phase_bwd)/(4*pi*dt);

% Compute first point using forward second-order formula
phase_01=angle(signal_col(2)*conj(signal_col(1)));
phase_12=angle(signal_col(3)*conj(signal_col(2)));
freq(1)=(3*phase_01-phase_12)/(4*pi*dt);

% Compute last point using backward second-order formula
npts=numel(signal_col);
phase_n1=angle(signal_col(npts)*conj(signal_col(npts-1)));
phase_n2=angle(signal_col(npts-1)*conj(signal_col(npts-2)));
freq(end)=(3*phase_n1-phase_n2)/(4*pi*dt);

% Restore the shape of the input vector
freq=reshape(freq,size(signal));

end

% Consistency enforcement
function grumble(signal,dt)
if (~isnumeric(signal))||(~isvector(signal))||isempty(signal)||isreal(signal)
    error('signal must be a non-empty complex numeric vector.');
end
if numel(signal)<3
    error('signal must contain at least three points.');
end
if any(~isfinite(signal))
    error('signal must not contain Inf or NaN.');
end
if (~isnumeric(dt))||(~isreal(dt))||(~isscalar(dt))||(~isfinite(dt))||(dt<=0)
    error('dt must be a positive real scalar.');
end
end


