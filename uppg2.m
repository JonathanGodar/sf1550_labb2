%% Interpolation

h = 0.25; % Får ej ändras i koden nedan

% Er kod här...
[t,x,y,vx,vy] = kastbana(h);

%% Linjär interpolation
x=x';
y=y';
degree = 1;

coefficient_matrix = piecewise_interpolation(x, y, degree);

plot_interpolation(x, coefficient_matrix, degree);

%%
calculate_and_plot_x_for_impact(y, coefficient_matrix, degree)

%% Kvadratisk interpolation

% Er kod här...


%%
lower_x = 1;
upper_x = 5;
grad = 1;

x_points = linspace(lower_x, upper_x, 123);
y_points = exp(x_points);

x_exact = linspace(lower_x, upper_x, 1000);
y_exact = exp(x_exact);

plot(x_exact, y_exact+10)
hold on

coeff = piecewise_interpolation(x_points, y_points, grad);
plot_interpolation(x_points, coeff, grad)




%%

lower_x = 1;
upper_x = 5;
grad = 1;

x_points = linspace(lower_x, upper_x, 21);
y_points = exp(x_points);

x_exact = linspace(lower_x, upper_x, 1000);
y_exact = exp(x_exact);

plot(x_exact, y_exact+10)
hold on


coeff = piecewise_interpolation(x_points, y_points, grad);
plot_interpolation(x_points, coeff, grad)
%%

% c_1, c_2 * x, c_3 *x^2 ...
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

function [coefficent_matrix] = piecewise_interpolation(x_points, y_points, grad)
    coefficent_matrix = [];
    for index = 1:(size(x_points,2)-1)/(grad)
        group_start = (index-1) * grad + 1;
        group_end = index*grad + 1;
        x_points_for_piece = x_points(group_start : group_end);
        y_points_for_piece = y_points(group_start : group_end);
        coefficents_for_piece = interpolate(x_points_for_piece,y_points_for_piece,grad);
        coefficent_matrix = [coefficent_matrix coefficents_for_piece];
    end
end

function plot_interpolation(x_points, coefficent_matrix, grad)
    hold on

    for index = 1:(size(x_points,2)-1)/grad
        group_start = (index-1) * grad + 1;
        group_end = index*grad + 1;
        x = linspace(x_points(group_start), x_points(group_end), 100);
        plot(x,evaluate_polynomial_at(coefficent_matrix(:,index)',x))
    end
end

function calculate_and_plot_x_for_impact(y_points, coefficent_matrix, grad)
    for index = (1:size(y_points,2)-1)
        if y_points(index) >= 0 && y_points(index+1) <=0
            indexs_for_impact = index;
            plot([0,50],[y_points(index),y_points(index)])
            plot([0,50],[y_points(index+1),y_points(index+1)])
        end
    end

    coefficents_for_relevant_interval = coefficent_matrix(:,index-1);
    if grad == 1
        disp("Bannan")
        c = coefficents_for_relevant_interval;
        x_for_impact = -c(1)/c(2);
    elseif grad == 2
        c = coefficents_for_relevant_interval;
        x_for_impact = (-c(2)+sqrt(c(2)^2-4*c(3)*c(1))/(x*c(3)))
    else
        error("Not available grade at the moment, wait for version 31.2")
    end

    plot(x_for_impact,0,'ro')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t,x,y,vx,vy]=kastbana(h)

%KASTBANA(H) beräknar banan för ett kast med en liten boll.
%
%   Dynamiken ges av en ODE som inkluderar effekten av luftmotståndet,
%
%      r'' = -g*ez-sigma*r'*|r'|/m.
%
%   Funktionen beräknar bollens position och hastighet vid
%   tidpunkter separerade med en given steglängd. Bollen kastas från
%   (X,Y)=(0,2) med hastigheten 30 m/s i 45 graders vinkel uppåt.
%
%   Syntax:
%
%   [T,X,Y,VX,VY] = KASTBANA(H)
%
%   H       - Steglängden mellan tidpunkterna.
%   T       - Vektor med tidpunkter där bollens position och hastighet beräknats.
%   X, Y    - Vektorer med bollens x- och y-koordinater vid tidpunkterna.
%   VX, VY  - Vektorer med bollens hastigheter i x- och y-led vid tidpunkterna.

%% Tennisboll, specifikationer

m = 56e-3;     % Massan (kg) = 56 gram
ra = 6.6e-2/2; % 6.6 cm in diameter

g=9.81;      % Tyngdaccelerationen (m/s^2)

rho=1.2;     % Luftens densitet (kg/m^3)
A=ra^2*pi;   % Kroppens tvärsnittsarea (m^2)
Cd=0.47;     % Luftmotståndskoefficient,
% "drag coefficient" (dimensionslös)
% Läs mer på http://en.wikipedia.org/wiki/Drag_coefficient

sigma = rho*A*Cd/2; % Totala luftmotståndet

T  = 5;      % Sluttid
v0 = 30;     % Utkasthastighet
al = pi/4;   % Utkastvinkel

% Begynnelsevärden

r0 = [0 2]';                   % Position
r1 = [v0*cos(al) v0*sin(al)]'; % Hastighet

% ODEns högerled

f = @(u) [u(3:4); -u(3:4)*norm(u(3:4),2)*sigma/m - [0;g]];  % RHS

u = [r0;r1];
U = u';
t = 0:h:T;

% Runge-Kutta 4

for tn=t(1:end-1)
    s1 = f(u);
    s2 = f(u + h/2*s1);
    s3 = f(u + h/2*s2);
    s4 = f(u + h*s3);
    u = u + h/6*(s1 + 2*s2 + 2*s3 + s4);
    U = [U; u'];
end

x  = U(:,1);
y  = U(:,2);
vx = U(:,3);
vy = U(:,4);

end

