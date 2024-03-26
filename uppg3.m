
load dollarkurs.mat
U = USDSEK';
N = length(U);
tt=(1:N)';
t = 1:1:N;

%% 3a Linjär modell
% Er kod här...

f = @(x) [1 x]
plot_and_check_model(t, U, f)


%% 3b Linjär + periodisk modell
% coefficients = least_square_custom(t, X, 1);
% Perioden ser ut att vara ca 430
L = 430;
f = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x /L)];
X = [plot_and_check_model(t, U, f); L]


%% 3c Icke-linjär modell
% Er kod här...

Xprev = 10 + X 
tau = 1e-8

while norm(Xprev-X) > tau
    Xprev = X 
    X = X - (J(X,U)'*J(X,U))^-1*J(X,U)'*F(X,U);
end

f = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x /L)];

X
model = make_model(X,U);


error = model - U;
grid on
plot(t, U);
hold on
plot(t, model);
plot(t, error);


calculate_average_square_error(model,U)
%%

function F_vec = F(X, y)
    F_vec = []
    for t = 1:size(y,2)
        row = [X(1) + X(2)*t + X(3) * sin(2*pi*t/X(5)) + X(4)*cos(2*pi*t/X(5)) - y(t)];
        F_vec = [F_vec; row];
    end
end

function F_vec = make_model(X, y)
    F_vec = []
    for t = 1:size(y,2)
        row = [X(1) + X(2)*t + X(3) * sin(2*pi*t/X(5)) + X(4)*cos(2*pi*t/X(5))];
        F_vec = [F_vec; row];
    end
end


function F_jacobian= J(X, y)
    F_jacobian = []
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
    
    calculate_average_square_error(model,y_values)
end


function [coeffs] = least_square(x_points, y_points, f) 
    val_matrix = [];
    for x_point = x_points  
        row =  f(x_point) % [1 x_point sin(2* pi * x_point / L) cos(2*pi * x_point /L)]
        val_matrix = [val_matrix; row];
    end
    coeffs = val_matrix\y_points';
end

% function [coeffs] = least_square_custom(x_points, y_points, grad) 
%     L = 430;
%     val_matrix = [];
%     for x_point = x_points  
%         row = [];
%         row = [1 x_point sin(2* pi * x_point / L) cos(2*pi * x_point /L)]
%         val_matrix = [val_matrix; row];
%     end
%     coeffs = val_matrix\y_points';
% end

function [coeffs] = least_square_polynomial(x_points, y_points, grad) 
    val_matrix = [];
    for x_point = x_points  
        row = [];
        for exponent = 0:grad
            row = [row x_point^exponent];
        end
        val_matrix = [val_matrix; row];
    end
    coeffs = val_matrix\y_points';
end

function evaluated = evaluate_polynomial_at(coefficients, x_values) 
    evaluated = [];
    for x_value = x_values
        sum = 0;
        for coeff_idx= [1:size(coefficients,2)]
            sum = sum + coefficients(coeff_idx) * x_value^(coeff_idx-1);
        end
        evaluated = [evaluated; sum];
    end
end

function average_square_error = calculate_average_square_error(actual_value,model_value)
    N = length(actual_value);
    sum_of_squares = 0;
    for index = (1:N)
        sum_of_squares = sum_of_squares + (actual_value(index) - model_value(index))^2;
    end

    average_square_error = sum_of_squares/N;
end