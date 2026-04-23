% Detailed unit test for ttclass/amensolve against dense references.
%
% ilya.kuprov@weizmann.ac.il
% d.savostyanov@soton.ac.uk
%
% The test builds structured positive-definite tensor-train linear systems,
% solves them with AMEn, and compares the result against dense direct solves
% and dense residuals. The cases include exact small systems, a dense-reference
% case of dimension 2000, a zero-enrichment regression, and a nonsymmetric
% finite-output smoke.

function amensolve_test_1()

% Initialise the random number generator
rng(1);

% Build the main accuracy test cases
cases={...
    struct('name','small_exact','dims',[5 4 3 2],'nterms_A',3,'nterms_y',2,'diag_shift',5e-2,'tol',1e-10,...
           'opts',struct('nswp',80,'init_guess_rank',2,'enrichment_rank',4,'rmax',24,'max_full_size',600,'local_iters',150,'verb',0),...
           'err_factor',200,'res_factor',50),...
    struct('name','medium_balanced','dims',[10 10 10],'nterms_A',3,'nterms_y',2,'diag_shift',8e-2,'tol',1e-8,...
           'opts',struct('nswp',120,'init_guess_rank',2,'enrichment_rank',6,'rmax',32,'max_full_size',500,'local_iters',220,'verb',0),...
           'err_factor',500,'res_factor',80),...
    struct('name','large_2000','dims',[20 10 10],'nterms_A',2,'nterms_y',2,'diag_shift',1e-1,'tol',2e-8,...
           'opts',struct('nswp',140,'init_guess_rank',2,'enrichment_rank',6,'rmax',24,'max_full_size',450,'local_iters',260,'verb',0),...
           'err_factor',600,'res_factor',100)};

% Run the main dense-reference tests
for n=1:numel(cases)

    % Pull out the current case
    case_data=cases{n};

    % Build the operator, right-hand side, and dense references
    [A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_spd_case(case_data.dims,case_data.nterms_A,case_data.nterms_y,case_data.diag_shift,n);

    % Run the AMEn solve
    x_tt=amensolve(A_tt,y_tt,case_data.tol,case_data.opts);
    x_amen=full(x_tt);

    % Compute the achieved relative solution error
    rel_err=norm(x_amen-x_dense,2)/max(norm_x,eps);

    % Check the residual against the dense system
    rel_res=norm(A_dense*x_amen-y_dense,2)/max(norm(y_dense,2),eps);

    % Practical acceptance thresholds for an iterative approximate solver
    err_lim=max(case_data.err_factor*case_data.tol,1e-10);
    res_lim=max(case_data.res_factor*case_data.tol,1e-10);

    % Check the dense-reference agreement
    if rel_err>err_lim
        error('amensolve_test_1(%s): relative solution error %.3e exceeds %.3e.',...
              case_data.name,rel_err,err_lim);
    end

    % Check the residual contract
    if rel_res>res_lim
        error('amensolve_test_1(%s): relative residual %.3e exceeds %.3e.',...
              case_data.name,rel_res,res_lim);
    end

    % Check that the result is a single tensor train
    if x_tt.ntrains~=1
        error('amensolve_test_1(%s): amensolve output is not a single tensor train.',...
              case_data.name);
    end

    % Check that the physical dimensions are preserved
    if any(x_tt.sizes(:,2)~=1) || any(x_tt.sizes(:,1)~=A_tt.sizes(:,2))
        error('amensolve_test_1(%s): output mode sizes do not match the linear system.',...
              case_data.name);
    end

    % Check that no non-finite numbers appeared
    if ~all(isfinite(x_amen(:)))
        error('amensolve_test_1(%s): non-finite values detected in the result.',...
              case_data.name);
    end

    % Report progress
    fprintf('amensolve_test_1: %s passed, dim=%d, relerr=%3.3e, relres=%3.3e, maxrank=%d\n',...
            case_data.name,numel(x_dense),rel_err,rel_res,max(x_tt.ranks));

end

% Reproducibility check on a medium problem
[A_rep,y_rep,~,~,x_rep,norm_rep]=build_spd_case([10 10 10],3,2,8e-2,99);
opts_rep=struct('nswp',120,'init_guess_rank',2,'enrichment_rank',6,'rmax',32,'max_full_size',500,'local_iters',220,'verb',0);
rng(12345);
x_one=amensolve(A_rep,y_rep,1e-8,opts_rep);
rng(12345);
x_two=amensolve(A_rep,y_rep,1e-8,opts_rep);
err_one=norm(full(x_one)-x_rep,2)/max(norm_rep,eps);
err_two=norm(full(x_two)-x_rep,2)/max(norm_rep,eps);
if abs(err_one-err_two)>1e-11
    error('amensolve_test_1: reproducibility check failed, relative errors differ.');
end
fprintf('amensolve_test_1: reproducibility passed, relerr=%3.3e\n',err_one);

% Zero-enrichment regression on a dense-reference system
[A_zero,y_zero,A_zero_dense,y_zero_dense,x_zero_dense,norm_zero]=build_spd_case([8 10 10],3,2,8e-2,77);
x_zero=amensolve(A_zero,y_zero,1e-8,struct('nswp',140,'init_guess_rank',2,'enrichment_rank',0,'rmax',32,'max_full_size',400,'local_iters',220,'verb',0));
x_zero_full=full(x_zero);
err_zero=norm(x_zero_full-x_zero_dense,2)/max(norm_zero,eps);
res_zero=norm(A_zero_dense*x_zero_full-y_zero_dense,2)/max(norm(y_zero_dense,2),eps);
if (~isfinite(err_zero))||(~isfinite(res_zero))
    error('amensolve_test_1: zero-enrichment regression produced non-finite output.');
end
fprintf('amensolve_test_1: zero_enrichment passed, dim=%d, relerr=%3.3e, relres=%3.3e, maxrank=%d\n',...
        numel(x_zero_dense),err_zero,res_zero,max(x_zero.ranks));

% Nonsymmetric smoke test: this path is not the main contract of amensolve,
% but it should still return a finite answer on a modest dense-reference case.
[A_ns_tt,y_ns_tt,A_ns_dense,y_ns_dense,~,~]=build_nonsym_case([6 6 5],4,4,123);
x_ns=amensolve(A_ns_tt,y_ns_tt,1e-8,struct('nswp',80,'init_guess_rank',2,'enrichment_rank',2,'rmax',32,'max_full_size',300,'local_iters',120,'verb',0));
x_ns_full=full(x_ns);
res_ns=norm(A_ns_dense*x_ns_full-y_ns_dense,2)/max(norm(y_ns_dense,2),eps);
if (~all(isfinite(x_ns_full(:))))||(~isfinite(res_ns))
    error('amensolve_test_1: nonsymmetric smoke test produced a non-finite result.');
end
fprintf('amensolve_test_1: nonsymmetric smoke passed, dim=%d, relres=%3.3e, maxrank=%d\n',...
        numel(x_ns_full),res_ns,max(x_ns.ranks));

% Final report
fprintf('amensolve_test_1: all tests passed.\n');

end

function [A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_spd_case(dims,nterms_A,nterms_y,diag_shift,seed)

% Reset the random number generator for this case
rng(seed);

% Build a buffered operator factor B and a buffered exact solution x
B_terms=build_rectangular_terms(dims,nterms_A);
x_terms=build_vector_terms(dims,nterms_y);

% Instantiate tensor trains
B_tt=ttclass(ones(1,nterms_A),B_terms,zeros(1,nterms_A));
x_exact_tt=ttclass(ones(1,nterms_y),x_terms,zeros(1,nterms_y));

% Form a Hermitian positive-definite system and the right-hand side
A_tt=shrink(B_tt'*B_tt+diag_shift*unit_like(B_tt'*B_tt));
y_tt=shrink(A_tt*x_exact_tt);

% Dense references
A_dense=full(A_tt);
x_dense=full(x_exact_tt);
y_dense=A_dense*x_dense;
norm_x=norm(x_dense,2);

end

function [A_tt,y_tt,A_dense,y_dense,x_dense,norm_x]=build_nonsym_case(dims,nterms_A,nterms_y,seed)

% Reset the random number generator for this case
rng(seed);

% Build operator and exact solution factors
A_terms=build_rectangular_terms(dims,nterms_A);
x_terms=build_vector_terms(dims,nterms_y);

% Instantiate tensor trains
A_base=ttclass(ones(1,nterms_A),A_terms,zeros(1,nterms_A));
x_exact_tt=ttclass(ones(1,nterms_y),x_terms,zeros(1,nterms_y));

% Add a small identity component to keep the matrix well-conditioned enough
A_tt=shrink(A_base+0.5*unit_like(A_base));
y_tt=shrink(A_tt*x_exact_tt);

% Dense references
A_dense=full(A_tt);
x_dense=full(x_exact_tt);
y_dense=A_dense*x_dense;
norm_x=norm(x_dense,2);

end

function terms=build_rectangular_terms(dims,nterms)

% Read the number of cores
ncores=numel(dims);

% Preallocate storage
terms=cell(ncores,nterms);

% Build all Kronecker factors
for n=1:nterms
    for k=1:ncores

        % Random complex matrix with controlled spectral scale
        current_factor=randn(dims(k))+1i*randn(dims(k));
        current_factor=current_factor/sqrt(dims(k));
        if n==k
            current_factor=current_factor+eye(dims(k));
        end
        terms{k,n}=current_factor;

    end
end

end

function terms=build_vector_terms(dims,nterms)

% Read the number of cores
ncores=numel(dims);

% Preallocate storage
terms=cell(ncores,nterms);

% Build all Kronecker factors
for n=1:nterms
    for k=1:ncores

        % Random complex vector factor with controlled norm
        current_factor=randn(dims(k),1)+1i*randn(dims(k),1);
        current_factor=current_factor/max(norm(current_factor,2),eps);
        terms{k,n}=current_factor;

    end
end

end
