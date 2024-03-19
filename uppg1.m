%% 1b -- Visualisering

load eiffel1


[eigen_values_matrix, eigen_vectors_matrix] = eig(A);
eigen_values = diag(eigen_vectors_matrix);

[eigen_values, sort_order] = sort(eigen_values);

for i = 1:4 % vector_idx = sort_order[1,4]
	figure(i);
	vector_idx = sort_order(i);
	y = eigen_vectors_matrix(:, vector_idx) * 1;
	trussplot(xnod + y(1:2:end), ynod + y(2:2:end), bars);
	break
end

%% Animering

y = eigen_vectors_matrix(:, sort_order(2))
trussanim(xnod, ynod, bars, y);

%% 1c -- Ber�kning av st�rsta och minsta egenv�rdena

% Er kod h�r...


%% 1d -- Ber�kning av andra egenv�rden

% Er kod h�r...


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = potens(A,tau)
%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (sk�l�r)
%
%  Utdata:
%
%  mu - st�rsta egenv�rdet till A (skal�r)
%  iter - antal iterationer som anv�nts (skal�r)

	y = rand(size(A, 2), 1)
	iter = 0
	mu = 0
	mu_prev = tolerance * 20
	while abs(mu - mu_prev) > tau
		iter = iter + 1;
		v = A * y;
		mu_prev = mu;
		mu = y' * v;
		y = v / norm(v);
	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mu, iter] = inverspotens(A,tau)

%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (sk�l�r)
%
%  Utdata:
%
%  mu - minsta egenv�rdet till A (skal�r)
%  iter - antal iterationer som anv�nts (skal�r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Er kod h�r...

end
