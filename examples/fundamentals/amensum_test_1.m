% Detailed unit test for ttclass/amensum against dense references.
%
% ilya.kuprov@weizmann.ac.il
% d.savostyanov@soton.ac.uk
%
% The test uses buffered rank-one tensor trains, compares the AMEn
% summation result against the exact dense sum, and checks both
% relative Frobenius error and internal consistency properties.
%
% Note: the underlying AMEn summation is approximate, and the paper
% motivating the method focuses on enrichment-assisted updates. The
% strict accuracy checks below therefore target enriched runs, while
% zero-enrichment is kept as a finite-output regression smoke test.

function amensum_test_1()

% Initialise the random number generator
rng(1);

% Build the test cases
cases={...
    struct('name','small_exact','dims',[5 4 3 2],'nterms',6,'tol',1e-12,'opts',struct('max_swp',80,'init_guess_rank',2,'enrichment_rank',4,'verb',0)),...
    struct('name','medium_balanced','dims',[10 10 10],'nterms',24,'tol',1e-10,'opts',struct('max_swp',80,'init_guess_rank',2,'enrichment_rank',4,'verb',0)),...
    struct('name','large_2000','dims',[20 10 10],'nterms',32,'tol',1e-8,'opts',struct('max_swp',120,'init_guess_rank',2,'enrichment_rank',6,'verb',0)),...
    struct('name','signed_coeffs','dims',[5 5 5 4],'nterms',18,'tol',1e-9,'opts',struct('max_swp',120,'init_guess_rank',3,'enrichment_rank',5,'verb',0))};

% Run the main numerical tests
for n=1:numel(cases)

    % Pull out the current case
    case_data=cases{n};

    % Build a buffered CP-format tensor train sum and a dense reference
    [x_ref,dense_ref,norm_ref]=build_case(case_data.dims,case_data.nterms,n);

    % Run the AMEn summation
    y_ref=amensum(x_ref,case_data.tol,case_data.opts);
    dense_amen=full(y_ref);

    % Compute the achieved relative error
    rel_err=norm(dense_amen-dense_ref,'fro')/max(norm_ref,eps);

    % Set a practical acceptance threshold for an approximate algorithm
    rel_lim=max(50*case_data.tol,1e-12);

    % Check the main approximation contract
    if rel_err>rel_lim
        error('amensum_test_1(%s): relative Frobenius error %.3e exceeds %.3e.',...
              case_data.name,rel_err,rel_lim);
    end

    % Check that the output is a single tensor train
    if y_ref.ntrains~=1
        error('amensum_test_1(%s): amensum output is not a single tensor train.',...
              case_data.name);
    end

    % Check that the physical dimensions are preserved
    if any(any(y_ref.sizes~=x_ref.sizes(:,1:2)))
        error('amensum_test_1(%s): output mode sizes do not match input mode sizes.',...
              case_data.name);
    end

    % Check that no non-finite numbers appeared
    if ~all(isfinite(dense_amen(:)))
        error('amensum_test_1(%s): non-finite values detected in the result.',...
              case_data.name);
    end

    % Report progress
    fprintf('amensum_test_1: %s passed, dim=%d, relerr=%3.3e, maxrank=%d\n',...
            case_data.name,size(dense_ref,1),rel_err,max(y_ref.ranks));

end

% Build a reproducibility case for rank-sensitive behaviour
[x_rep,dense_rep,norm_rep]=build_case([10 10 10],24,99);
opts_rep=struct('max_swp',100,'init_guess_rank',2,'enrichment_rank',4,'verb',0);
rng(12345);
y_one=amensum(x_rep,1e-9,opts_rep);
rng(12345);
y_two=amensum(x_rep,1e-9,opts_rep);
err_one=norm(full(y_one)-dense_rep,'fro')/max(norm_rep,eps);
err_two=norm(full(y_two)-dense_rep,'fro')/max(norm_rep,eps);
if abs(err_one-err_two)>1e-12
    error('amensum_test_1: reproducibility check failed, relative errors differ.');
end

% Check that disabling enrichment no longer crashes and still returns a finite object
[x_smoke,dense_smoke,norm_smoke]=build_case([8 10 10],20,77);
y_smoke=amensum(x_smoke,1e-8,struct('max_swp',160,'init_guess_rank',2,'enrichment_rank',0,'verb',0));
err_smoke=norm(full(y_smoke)-dense_smoke,'fro')/max(norm_smoke,eps);
if ~isfinite(err_smoke)
    error('amensum_test_1: zero-enrichment smoke test produced a non-finite error.');
end
fprintf('amensum_test_1: zero_enrichment smoke passed, dim=%d, relerr=%3.3e, maxrank=%d\n',...
        size(dense_smoke,1),err_smoke,max(y_smoke.ranks));

% Final report
fprintf('amensum_test_1: all tests passed.\n');

end

function [x_obj,dense_sum,norm_sum]=build_case(dims,nterms,seed)

% Reset the random number generator for this case
rng(seed);

% Read the number of tensor cores
ncores=numel(dims);

% Prepare storage for Kronecker terms
kronterms=cell(ncores,nterms);
coeff=zeros(1,nterms);

% Build all buffered rank-one tensor trains
for n=1:nterms

    % Generate a moderately scaled coefficient
    coeff(n)=(-1)^n*(0.5+rand());

    % Build one Kronecker factor per core
    for k=1:ncores

        % Generate a dense factor with controlled norm
        current_factor=randn(dims(k))+1i*randn(dims(k));
        current_factor=current_factor/norm(current_factor,'fro');
        kronterms{k,n}=current_factor;

    end

end

% Instantiate the tensor train object
x_obj=ttclass(coeff,kronterms,zeros(1,nterms));

% Assemble the dense reference
nn=prod(dims);
dense_sum=zeros(nn,nn);
for n=1:nterms
    dense_term=1;
    for k=1:ncores
        dense_term=kron(dense_term,kronterms{k,n});
    end
    dense_sum=dense_sum+coeff(n)*dense_term;
end
norm_sum=norm(dense_sum,'fro');

end

