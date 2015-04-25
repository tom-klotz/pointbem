function b = besselBackwardRecurrence(bn, bnplus1, M,z)

b = zeros(1 + M,1);
b(1+ M) = bnplus1;
b(1+ M-1) = bn;
for n=M-1:-1:1
    b(1+ (n-1)) = ((2*n+1)/z) * b(1+ n) - b(1+ (n+1));
end
