Iexact = 6.231467927023725;  % Ett noggrannt värde för I

%% 5a Trapetsregeln i 10 dimensioner

% Er kod här...
n = 8;
tic;
I = trapets10d(n);
time = toc;
error = abs(Iexact - I);
disp(["Approxiamtivt I med n = " + n + ": " + I; "Felet: " + error ; "Beräkningstid: " + time]);
%"Approxiamtivt I med n = 7: 6.2318"
% "Felet: 0.00029003"
% "Beräkningstid: 9.8023"

% "Approxiamtivt I med n = 8: 6.2317"
% "Felet: 0.00021863"
% "Beräkningstid: 34.4759"



%% 5b Monte-Carlo
N = 1e7;

x = 1:N;

errors = [];
list_of_list_of_approximated_integrals = [];
time = [];
for i = 1:5
    tic
    approximated_integrals = montecarlo(N);
    time = [time toc];
    list_of_approximated_integrals = approximated_integrals(1:N);
    list_of_list_of_approximated_integrals = [list_of_list_of_approximated_integrals list_of_approximated_integrals'];

    error = abs(Iexact - list_of_approximated_integrals)';
    errors = [errors, error];
    
end

average_error = sum(errors,2)/5;

figure(1);
plot(x, list_of_list_of_approximated_integrals);
axis([20 N 6.22 6.24])

figure(2);
loglog(x, errors);
hold on
loglog(x, x.^-.5 / x(end)^-.5 * average_error(end), "--", 'LineWidth', 8);
axis([20 N 1e-5 1])
% loglog(x, x.^-.5);
figure(3)
loglog(x,average_error);
hold on
loglog(x, x.^-.5 / x(end)^-.5 * average_error(end));
axis([20 N 1e-5 1])


format long
disp("Bästa gissning: ");
disp(list_of_approximated_integrals(end));
disp("Väntat fel:");
disp(1/sqrt(N));
disp("Faktisk fel:");
disp(error(end));
disp("Tid: " + time)



%%
tic;
result = montecarlo(1e6 + 5e5);
disp(toc);

disp(result(end) - Iexact)

%%

% Er kod här...
%
function I = montecarlo(N_max)
    f = @(x) exp(prod(x));
    points = rand(10, N_max) * 1.2;
    size(points)

    size_omega = 1.2^10;
    
    running_sum = zeros(1,N_max);
    running_sum(1) = f(points(1));
    i = 2;
    for point = points(:, 2:end)
        running_sum(i) = running_sum(i-1) + f(point);
        % running_sum = [running_sum, running_sum(end) + f(point)];
        i = i + 1;
    end
    
    n = 1:N_max;
    disp("Run " + size(running_sum))
    disp("sise " + size(n))
    I = size_omega ./ n .* running_sum;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I = trapets10d(n)
%  Indata:
%
%  n  - antal delintervall i varje koordinatriktning (skalär)
%
%  Utdata:
%
%  I - integralvärdet (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=1.2;

h = L/n;

x = 0:h:L;

I1 = zeros(n+1,1);
I2 = zeros(n+1,1);
I3 = zeros(n+1,1);
I4 = zeros(n+1,1);
I5 = zeros(n+1,1);
I6 = zeros(n+1,1);
I7 = zeros(n+1,1);
I8 = zeros(n+1,1);
I9 = zeros(n+1,1);
I10 = zeros(n+1,1);

for j1=0:n
    for j2=0:n
        for j3=0:n
            for j4=0:n
                for j5=0:n
                    for j6=0:n
                        for j7=0:n
                            for j8=0:n
                                for j9=0:n
                                    for j10=0:n
                                        I10(j10+1) = exp(j1*j2*j3*j4*j5*j6*j7*j8*j9*j10*h^10);
                                    end
                                    I9(j9+1) = h*(sum(I10) - I10(1)/2 - I10(end)/2);
                                end
                                I8(j8+1) = h*(sum(I9) - I9(1)/2 - I9(end)/2);
                            end
                            I7(j7+1) = h*(sum(I8) - I8(1)/2 - I8(end)/2);
                        end
                        I6(j6+1) = h*(sum(I7) - I7(1)/2 - I7(end)/2);
                    end
                    I5(j5+1) = h*(sum(I6) - I6(1)/2 - I6(end)/2);
                end
                I4(j4+1) = h*(sum(I5) - I5(1)/2 - I5(end)/2);
            end
            I3(j3+1) = h*(sum(I4) - I4(1)/2 - I4(end)/2);
        end
        I2(j2+1) = h*(sum(I3) - I3(1)/2 - I3(end)/2);
    end
    I1(j1+1) = h*(sum(I2) - I2(1)/2 - I2(end)/2);
end
I = h*(sum(I1) - I1(1)/2 - I1(end)/2);

end
