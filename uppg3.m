
load dollarkurs.mat
X = USDSEK';
N = length(X);
tt=(1:N)';

%% 3a Linjär modell

% Er kod här...
day = 1:1:N;
plot(day, USDSEK)
hold on
coefficients = interpolate(day, X, 1);

models = evaluate_polynomial_at(coefficients', day)

plot(day, X)
plot(day, )
error = 
plot(day, evaluate_polynomial_at(coefficients', day))

%% 3b Linjär + periodisk modell

% Er kod här...


%% 3c Icke-linjär modell

% Er kod här...

function [coeffs] = interpolate(x_points, y_points, grad) 
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

function average_square_error = calculate_average_square_error (t_value,function_value)
    N = length(t_value);
    sum_of_squares = 0;
    for index = (1:N)
        sum_of_squares = sum_of_squares + (t_value(index) - function_value(index))^2
    end

    average_square_error = sum_of_squares/N
end