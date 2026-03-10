% Performs a one-step BFGS Hessian update for maximisation using
% the argument and gradient increments from the previous step.
% Syntax:
%
%                     H=bfgs_upd(H,dx,dg)
%
% Arguments:
%
%    H      - current BFGS approximation to the Hessian
%             matrix corresponding to the *negative*
%             Hessian of the objective; use [] on the
%             first call
%
%    dx     - increment in arguments between the current
%             and the previous step
%
%    dg     - increment in gradients between the current
%             and the previous step
%
% Returns:
%
%    H      - updated BFGS approximation to the Hessian
%             matrix corresponding to the *negative*
%             Hessian of the objective
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bfgs_upd.m>

function H=bfgs_upd(H,dx,dg)

% Check consistency
grumble(H,dx,dg);

% Convert increments into column vectors
dx=dx(:); dg=dg(:);

% Curvature tolerance
curv_rel_tol=0.01;

% Relevant inner products
dgdx=dg'*dx;
dgdg=dg'*dg;
dxdx=dx'*dx;

% Curvature pair validation
pair_ok=isfinite(dgdx)&&isfinite(dgdg)&&isfinite(dxdx)&&...
        (dgdg>0)&&(dxdx>0)&&(dgdx<-curv_rel_tol*sqrt(dgdg*dxdx));

% Return identity if first pair is bad
if isempty(H)&&(~pair_ok)
    H=eye(numel(dx));
    return;
end

% Return unchanged Hessian if pair is bad
if (~isempty(H))&&(~pair_ok)
    H=real((H+H')/2);
    return;
end

% Convert gradient increment sign for maximisation
y=-dg;

% Initialise Hessian estimate when absent
if isempty(H)

    % Relevant inner products
    ydx=y'*dx;
    yyy=y'*y;

    % Unit matrix scaling
    if (~isfinite(ydx))||(~isfinite(yyy))||(ydx<=eps)||(yyy<=0)
        gamma=1;
    else
        gamma=yyy/ydx;
    end

    % Initial pseudo-Hessian guess
    H=gamma*eye(numel(dx));

end

% Build BFGS components
ydx=y'*dx;
Hdx=H*dx;
dxHdx=dx'*Hdx;

% Return unchanged Hessian if update is numerically unsafe
if (~isfinite(ydx))||(ydx<=eps)
    H=real((H+H')/2);
    return;
end
if (~isfinite(dxHdx))||(dxHdx<=eps)
    H=real((H+H')/2);
    return;
end

% Plain BFGS Hessian update
H=H-(Hdx*Hdx')/dxHdx+(y*y')/ydx;

% Numerical clean-up
H=real((H+H')/2);

end

% Consistency enforcement
function grumble(H,dx,dg)
if (~isnumeric(dx))||(~isreal(dx))||(~isvector(dx))||isempty(dx)
    error('dx must be a non-empty real numeric vector.');
end
if (~isnumeric(dg))||(~isreal(dg))||(~isvector(dg))||isempty(dg)
    error('dg must be a non-empty real numeric vector.');
end
if numel(dx)~=numel(dg)
    error('dx and dg must have the same number of elements.');
end
if ~isempty(H)
    if (~isnumeric(H))||(~isreal(H))||(size(H,1)~=size(H,2))
        error('H must be empty or a real square matrix.');
    end
    if size(H,1)~=numel(dx)
        error('H dimension must match the number of elements in dx.');
    end
end
end

% The easiest way to solve a problem is to deny it exists.
%
% Isaac Asimov


