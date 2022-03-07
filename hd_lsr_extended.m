function [A, B, error] = hd_lsr_extended(Y, s, normalize)
% Algorithm 1 from the paper
[n, m] = size(Y);

if s <= m
    A = Y(:, 1:s)./abs(Y(: ,1:s));
else
    A = randn(size(Y))+1i*randn(size(Y));
    A = A./abs(A);
end

B = A\Y; B = [B; conj(B)];
A = [A zeros(size(A))];

ss = min(s, m);

error = [norm(Y - A*B, 'fro')^2/norm(Y, 'fro')^2];
K = 15; % iterations
for k = 1:K
    R = Y - A*B;
    for j = 1:ss
        R = R + A(:, [j ss+j])*B([j ss+j], :);
        for i = 1:n
%             R(i, :) = R(i, :) + A(i, [j s+j])*B([j s+j], :);
            
            %lineR = R(i,:); line1B = B(j,:);
            a = sum(real(R(i,:)).*real(B(j,:)));
            b = sum(real(R(i,:)).*imag(B(j,:)));
            c = sum(imag(R(i,:)).*real(B(j,:)));
            d = sum(imag(R(i,:)).*imag(B(j,:)));
            
            if a*d < b*c
                A(i, j) = 0;
                A(i, ss+j) = a-d+1i*(c+b);
                A(i, ss+j) = A(i, ss+j)/abs(A(i, ss+j));
            else
                A(i, j) = a+d+1i*(c-b);
                A(i, j) = A(i, j)/abs(A(i, j));
                A(i, ss+j) = 0;
            end

%             R(i, :) = R(i, :) - A(i, [j s+j])*B([j s+j], :);
        end
        R = R - A(:, [j ss+j])*B([j ss+j], :);
    end
    AA = [real(A(:, 1:ss)+A(:, ss+(1:ss))) -imag(A(:, 1:ss)-A(:, ss+(1:ss))); imag(A(:, 1:ss)+A(:, ss+(1:ss))) real(A(:, 1:ss)-A(:, ss+(1:ss)))]\[real(Y); imag(Y)];
    B = [AA(1:ss, :) + 1i*AA(ss+(1:ss), :); AA(1:ss, :) - 1i*AA(ss+(1:ss), :)];

    if normalize
      B = B*sqrt(s)/norm(A*B,'fro');
    end

    error = [error norm(Y - A*B, 'fro')^2/norm(Y, 'fro')^2];
end
