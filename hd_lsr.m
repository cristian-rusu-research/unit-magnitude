function [A, B, error] = hd_lsr(Y, s, normalize)
[n, m] = size(Y);

if s <= m
    A = Y(:, 1:s)./abs(Y(: ,1:s));
else
    A = randn(size(Y))+1i*randn(size(Y));
    A = A./abs(A);
end
B = A\Y;

ss = min(s, m);

error = [norm(Y - A*B, 'fro')^2/norm(Y, 'fro')^2];
K = 15; % iterations
for k = 1:K
    R = Y - A*B;
    for j = 1:ss
        R = R + A(:, j)*B(j, :);
        
        for i = 1:n
            x = sum(conj(B(j,:)).*R(i,:));
            A(i,j) = x/abs(x);
        end
        
        R = R - A(:, j)*B(j, :);
    end
    B = A\Y;
    % error = [error norm(Y - A*B, 'fro')^2/norm(Y, 'fro')^2];
    
    if normalize
        B = B*sqrt(s)/norm(A*B,'fro');
    end

    error = [error norm(Y - A*B, 'fro')^2/norm(Y, 'fro')^2];
end
