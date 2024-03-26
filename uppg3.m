
load dollarkurs.mat
X = USDSEK';
N = length(X);
tt=(1:N)';
days = 1:1:N;

%% 3a Linjär modell

% Er kod här...

hold on
coefficients = least_square_polynomial(days, X, 1);


model = evaluate_polynomial_at(coefficients', days)';
error = model - X;

grid on
plot(days, X); 
plot(days, model);
plot(days, error);

average_square_error = calculate_average_square_error(X,model)

%% 3b Linjär + periodisk modell
% coefficients = least_square_custom(days, X, 1);
L = 430;
f = @(x) [1 x sin(2* pi * x / L) cos(2*pi * x /L)];
plot_and_check_model(f)

% Perioden ser ut att vara ca 430

% Er kod här...


%% 3c Icke-linjär modell


function plot_and_check_model(f)
    coeffs = least_square(days, X, f);

    model = [];
    for day = days 
        model = [model, f(day)*coeffs];
    end
    % model = f_excplict(days)

    error = model - X;

    grid on
    plot(days, X);
    hold on
    plot(days, model);
    plot(days, error);
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