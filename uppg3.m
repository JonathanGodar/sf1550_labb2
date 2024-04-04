load dollarkurs.mat;
U = USDSEK';
N = length(U);
tt=(1:N)';
t = 1:1:N;

%% 3a Linjär modell
% Er kod här...

figure(1)
f = @(x) [1 x];
X = plot_and_check_model(t, U, f);
disp(["c_0 = ";"c_1 = "] + X);

%% 3b Linjär + periodisk modell
% coefficients = least_square_custom(t, X, 1);
% Perioden ser ut att vara ca 430
figure(2)
L = 430;
f = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x /L)];
X = [plot_and_check_model(t, U, f); L];
disp(["d_0 = ";"d_1 = ";"d_2 = ";"d_3 = "; "L = "] + X);


%% 3c Icke-linjär modell
% Er kod här...

figure(3)
Xprev = 10 + X;
tau = 1e-8;

while norm(Xprev-X) > tau
    Xprev = X;
    X = X - (J(X,U)'*J(X,U))^-1*J(X,U)'*F(X,U);
end

f = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x /L)];


model = make_model(X,U);

error = model' - U;
grid on
plot(t, U);

hold on
plot(t, model);

plot(t, error);

disp("Medelkvadratfelet = " + calculate_average_square_error(model,U));
disp(["d_0 = ";"d_1 = ";"d_2 = ";"d_3 = ";"L = "] + X);

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


function coeffs = plot_and_check_model(x_values, y_values, f)
    coeffs = least_square(x_values, y_values, f);

    model = [];
    for x = x_values 
        model = [model, f(x)*coeffs];
    end
    % model = f_excplict(y_values)
    error = model - y_values;

    grid on
    plot(x_values, y_values);
    hold on
    plot(x_values, model);
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