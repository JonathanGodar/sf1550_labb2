load dollarkurs.mat;
U = USDSEK';
N = length(U);
tt=(1:N)';
t = 1:1:N;

%% 3a Linjär modell
% Er kod här...

figure(1)
f_linear = @(x) [1 x];
coeffs_linear = plot_and_check_model(t, U, f_linear);
title("Linjär anpassning/Linjär anpassning fel");
disp(["c_0 = ";"c_1 = "] + coeffs_linear);

%% 3b Linjär + periodisk modell
% Perioden ser ut att vara ca 430
figure(2)
L = 430;
f_periodical = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x / L)];
coeffs_periodical = [plot_and_check_model(t, U, f_periodical); L];
title("Periodisk anpassning med gissning för L / Fel");
disp(["d_0 = ";"d_1 = ";"d_2 = ";"d_3 = "; "L = "] + coeffs_periodical);


%% 3c Icke-linjär modell
% Er kod här...

figure(3)
X = coeffs_periodical;
Xprev = 10 + X;
tau = 1e-8;

while norm(Xprev-X) > tau
    Xprev = X;
    X = X - (J(X,U)'*J(X,U))^-1*J(X,U)'*F(X,U);
end

f = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x /L)];


model = make_model(X,U);

error = model' - U;
subplot(2, 1, 1);
grid on
plot(t, U);

hold on
plot(t, model);

subplot(2, 1, 2);
plot(t, error);
title("Periodisk anpassning med beräkning av L/Fel")
disp("Medelkvadratfelet = " + calculate_average_square_error(model,U));
disp(["d_0 = ";"d_1 = ";"d_2 = ";"d_3 = ";"L = "] + X);

figure(4);
plot(t, U, t, model, t, evaluate_function_at(t, f_linear, coeffs_linear), t, evaluate_function_at(t, f_periodical, coeffs_periodical(1:end-1)));
legend("$", "Model 3", "Model 1", "Model 2");
%%


function F_vec = F(X, y)
    F_vec = [];
    for t = 1:size(y,2)
        row = [X(1) + X(2)*t + X(3) * sin(2*pi*t/X(5)) + X(4)*cos(2*pi*t/X(5)) - y(t)];
        F_vec = [F_vec; row];
    end
end

function F_vec = make_model(X, y)
    F_vec = [];
    for t = 1:size(y,2)
        row = [X(1) + X(2)*t + X(3) * sin(2*pi*t/X(5)) + X(4)*cos(2*pi*t/X(5))];
        F_vec = [F_vec; row];
    end
end


function F_jacobian= J(X, y)
    F_jacobian = [];
    for t = 1:size(y,2)
        row = [1 t  sin(2*pi*t/X(5)) cos(2*pi*t/X(5)) (X(4)*sin(2*pi*t/X(5)) - X(3)*cos(2*pi*t/X(5))) * 2*pi*t/(X(5)^2) ];
        F_jacobian = [F_jacobian; row];
    end
end


function y = evaluate_function_at(x_values, f, coeffs) 
    y = [];
    for x = x_values 
        y = [y, f(x) * coeffs];
    end
end

function coeffs = plot_and_check_model(x_values, y_values, f)
    coeffs = least_square(x_values, y_values, f);

    model = evaluate_function_at(x_values, f, coeffs);
    error = model - y_values;

    subplot(2, 1, 1);
    grid on
    plot(x_values, y_values);
    hold on
    plot(x_values, model);
    
    subplot(2, 1, 2);
    plot(x_values, error);
    
    MSE = calculate_average_square_error(model,y_values);
    disp("Medelkvadratfelet: " + MSE);
end


function [coeffs] = least_square(x_points, y_points, f) 
    val_matrix = [];
    for x_point = x_points  
        row =  f(x_point);
        val_matrix = [val_matrix; row];
    end
    coeffs = val_matrix\y_points';
end

function average_square_error = calculate_average_square_error(actual_value,model_value)
    N = length(actual_value);
    sum_of_squares = 0;
    for index = (1:N) 
        sum_of_squares = sum_of_squares + (actual_value(index) - model_value(index))^2;
    end

    average_square_error = sum_of_squares/N;
end