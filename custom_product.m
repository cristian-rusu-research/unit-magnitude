function RES = custom_product(Z, Y, MASK)
[n, p] = size(Z);
[~, m] = size(Y);

RES = zeros(n, m);
for i = 1:n
    for j = 1:m
        for k = 1:p
            if (MASK(i, k) == 0)
                RES(i, j) = RES(i, j) + conj(Z(i, k)*Y(k, j));
            else
                RES(i, j) = RES(i, j) + Z(i, k)*Y(k, j);
            end
        end
    end
end
