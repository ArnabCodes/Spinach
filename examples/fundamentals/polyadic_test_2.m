% Unit tests for advanced polyadic functionality.
%
% ilya.kuprov@weizmann.ac.il

function polyadic_test_2()

% Get random test matrices
a=randn(2,2)+1i*randn(2,2);
b=randn(3,3)+1i*randn(3,3);
c=sprandn(2,2,0.75)+1i*sprandn(2,2,0.75);
d=randn(3,3)+1i*randn(3,3);
left=randn(5,6)+1i*randn(5,6);
right=randn(6,4)+1i*randn(6,4);
pre=randn(8,6)+1i*randn(8,6);
suf=randn(6,7)+1i*randn(6,7);

% Build a reference polyadic and its matrix form
p=polyadic({{a,b},{c,d}});
ref=kron(a,b)+kron(c,d);

% Check constructor, full, inflate, and validate
assert_small(norm(full(p)-full(ref),1),1e-12);
assert_small(norm(full(inflate(p)-ref),1),1e-12);
validate(p);
disp('Constructor/full/inflate/validate test PASSED.');

% Check prefixes, suffixes, size, and emptiness
p_pref=prefix(left,p);
p_pref=suffix(p_pref,right);
ref_pref=left*ref*right;
assert_small(norm(full(p_pref)-full(ref_pref),1),1e-12);
assert_small(norm(size(p_pref)-size(ref_pref)),1e-12);
assert(isequal(isempty(p_pref),isempty(ref_pref)));
disp('Prefix/suffix/size/isempty test PASSED.');

% Check addition and subtraction paths
q=polyadic({{d.',a.'}});
ref_q=kron(d.',a.');
assert_small(norm(full(p+q)-full(ref+ref_q),1),1e-12);
assert_small(norm(full(p-q)-full(ref-ref_q),1),1e-12);
assert_small(norm(full(p+2)-full(ref+2),1),1e-12);
disp('Plus/minus test PASSED.');

% Check multiplication paths
assert_small(norm(full(3*p)-full(3*ref),1),1e-12);
assert_small(norm(full(p*3)-full(ref*3),1),1e-12);
assert_small(norm(full(pre*p)-full(pre*ref),1),1e-12);
assert_small(norm(full(p*suf)-full(ref*suf),1),1e-12);
assert_small(norm(full(sparse(pre)*p)-full(sparse(pre)*ref),1),1e-12);
assert_small(norm(full(p*sparse(suf))-full(ref*sparse(suf)),1),1e-12);
r=polyadic({{randn(2,4)+1i*randn(2,4),randn(3,2)+1i*randn(3,2)}});
ref_r=kron(r.cores{1}{1},r.cores{1}{2});
assert_small(norm(full(p*r)-full(ref*ref_r),1),1e-12);
disp('Mtimes test PASSED.');

% Check Kronecker products
k_num=randn(2,2)+1i*randn(2,2);
assert_small(norm(full(kron(p,k_num))-full(kron(ref,k_num)),1),1e-12);
assert_small(norm(full(kron(k_num,p))-full(kron(k_num,ref)),1),1e-12);
assert_small(norm(full(kron(p,q))-full(kron(ref,ref_q)),1),1e-12);
disp('Kron test PASSED.');

% Check transpose operations
assert_small(norm(full(transpose(p))-full(transpose(ref)),1),1e-12);
assert_small(norm(full(ctranspose(p))-full(ctranspose(ref)),1),1e-12);
disp('Transpose/ctranspose test PASSED.');

% Check finiteness and internal non-zero counts
assert(isequal(allfinite(p),all(isfinite(ref(:)))));
p_nan=prefix(NaN,p);
assert(~allfinite(p_nan));
assert(isequal(nnz(p),nnz(a)+nnz(b)+nnz(c)+nnz(d)));
disp('Allfinite/nnz test PASSED.');

% Check zero-dimension behaviour
p_empty=polyadic({{zeros(0,2)}});
assert(isempty(p_empty));
assert(isequal(size(p_empty),[0 2]));
disp('Empty-matrix test PASSED.');

% Check nested simplification paths
p_nested=polyadic({{polyadic({{eye(2)}}),b}});
s_nested=simplify(p_nested);
assert_small(norm(full(s_nested)-full(kron(eye(2),b)),1),1e-12);
p_nested_prefix=prefix(polyadic({{eye(6)}}),p);
s_nested_prefix=simplify(p_nested_prefix);
assert_small(norm(full(s_nested_prefix)-full(ref),1),1e-12);
p_nested_suffix=suffix(p,polyadic({{eye(6)}}));
s_nested_suffix=simplify(p_nested_suffix);
assert_small(norm(full(s_nested_suffix)-full(ref),1),1e-12);
disp('Simplify test PASSED.');

% Check GPU upload path when hardware is available
if gpuDeviceCount>0
    pg=gpuArray(p);
    assert_small(norm(gather(full(pg))-full(ref),1),1e-10);
    disp('GpuArray test PASSED.');
else
    disp('GpuArray test SKIPPED (no GPU available).');
end

end

function assert_small(value,threshold)
if value>=threshold
    error('numeric tolerance test failed.');
end
end


