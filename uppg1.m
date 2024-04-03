%% 1b -- Visualisering
format long
load eiffel1

[eigen_vectors_matrix, eigen_values_matrix] = eig(A);

eigen_values = diag(eigen_values_matrix);
[eigen_values, sort_order] = sort(eigen_values);

for i = 1:4 
	figure(i);
	vector_idx = sort_order(i);
	y = eigen_vectors_matrix(:, vector_idx)*3; 
	trussplot(xnod + y(1:2:end), ynod + y(2:2:end), bars);
end

%% Animering
y = eigen_vectors_matrix(:, sort_order(5));
trussanim(xnod, ynod, bars, y*1);
%% 1c -- Beräkning av största och minsta egenvärdena
tau = 1e-10;
load eiffel1;
[mini_mu_1, mini_iter_1] = potens(A, tau);
[maxi_mu_1, maxi_iter_1] = inverspotens(A, tau);


[eigen_values_matrix, eigen_vectors_matrix] = eig(A);
eigen_values = diag(eigen_vectors_matrix);
[eigen_values, sort_order] = sort(eigen_values);

load eiffel2;
[mini_mu_2, mini_iter_2] = potens(A, tau);
[maxi_mu_2, maxi_iter_2] = inverspotens(A, tau);

load eiffel3;
[mini_mu_3, mini_iter_3] = potens(A, tau);
[maxi_mu_3, maxi_iter_3] = inverspotens(A, tau);

load eiffel4;
[mini_mu_4, mini_iter_4] = potens(A, tau);
[maxi_mu_4, maxi_iter_4] = inverspotens(A, tau);

%% 
load eiffel1
disp("Eiffel 1");
eiffel_data_row = create_data_row(A, tau);

table = [1 eiffel_data_row];

load eiffel2
disp("Eiffel 2");
eiffel_data_row = create_data_row(A, tau);
table = [table; 2 eiffel_data_row];

load eiffel3
disp("Eiffel 3");
eiffel_data_row = create_data_row(A, tau);
table = [table; 3 eiffel_data_row];

load eiffel4
disp("Eiffel 4");
eiffel_data_row = create_data_row(A, tau);
table = [table; 4 eiffel_data_row];

disp(table);

%% 1d -- Beräkning av andra egenvärden

load eiffel1
tau = 1e-10;
sigmaArray = [10, 50, 67];
for sigma = sigmaArray
	[mu,iter] = inverspotens(A - sigma * eye(size(A,1)), tau);
	closestEigenvalueToSigma = mu + sigma;
	disp([sigma, closestEigenvalueToSigma, iter])
end

% disp(sort(eig(A)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [mu, iter] = potens(A,tau)
%  Indata:
%
%  A  - matrisen (kvadratisk)
%  tau - feltolerans (skalär)
%
%  Utdata:
%
%  mu - största egenvärdet till A (skalär)
%  iter - antal iterationer som använts (skalär)
	y =  rand(size(A, 2), 1); %ones(size(A, 2),1) * 4;
	A = sparse(A);
	iter = 0;
	mu = 0;
	mu_prev = tau * 20;

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
%  tau - feltolerans (skalär)
%
%  Utdata:
%
%  mu - minsta egenvärdet till A (skalär)
%  iter - antal iterationer som använts (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	A = inv(A);
	[mu, iter] = potens(A, tau);
	mu = 1/mu;
end


function data_row = create_data_row(A, tau)
	tic;
	[maxi_mu, maxi_iter] = potens(A, tau);
	[mini_mu, mini_iter] = inverspotens(A, tau);
	time_for_our_methods = toc;
	
	tic;
	eigen_values = eig(A);
	time_for_eig = toc;

	eigen_values = flip(sort(eigen_values));

	l2_over_l1 = eigen_values(2) / eigen_values(1);
	eig_size = size(eigen_values);
	eig_len = eig_size(1);
	ln_over_lnminus_one = eigen_values(eig_len)/eigen_values(eig_len-1);

	disp(["Eig time" "Our time"; time_for_eig time_for_our_methods]);
	disp(["Eig min/max" "Our min/max"; eigen_values(end) mini_mu; eigen_values(1) maxi_mu]);

	data_row = [maxi_mu maxi_iter l2_over_l1 mini_mu mini_iter ln_over_lnminus_one ];
end