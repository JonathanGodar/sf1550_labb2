

%% 4a Trapetsregeln och Simpson

format long
trap_30 = trapets(30)
trap_60 = trapets(60)

simp_30 = simpson(30)
simp_60 = simpson(60)

format long
disp(["intervall\Metod" "Trapets" "Simpson"; 
			"n = 30:" trap_30 simp_30; 
			"n = 60:" trap_60 simp_60]);

%% 4b Konvergensstudier

% Er kod här...
%1
ref = simpson(1e5);
err_trap = [];
x_max = 1000

x_trap = 1:x_max;
for n = x_trap
    err_trap = [err_trap abs(ref - trapets(n))];
end

x_simp = 2:2:x_max;
err_simp = [];
for n = x_simp
    err_simp = [err_simp abs(ref - simpson(n))];
end

R = 3;
loglog(R*x_trap.^-1,err_trap);
hold on;
loglog(R*x_simp.^-1,err_simp);
grid on

x_comp_trap = linspace(R/x_trap(1),R/x_trap(end),100);
y_comp_trap = x_comp_trap.^2 / x_comp_trap(1) * err_trap(1);
loglog(x_comp_trap, y_comp_trap)

x_comp_simp = linspace(R/x_simp(1),R/x_simp(end),100);
y_comp_simp = x_comp_simp.^4 / x_comp_simp(1) * err_simp(1);
loglog(x_comp_simp, y_comp_simp)
% y1 = x.^3 / sizeList(1)^3 * timeList(1); % +   % timeList(1);
% y2 = x.^2 / sizeList(1)^2 * timeList(1);
% y3 = x.^2.475 / sizeList(1)^2.475 * timeList(1);


%% 2
R=3;

series_of_n = 60:2:80
ratio_trap = [];
for n = series_of_n 
    ratio = (trapets(n)-trapets(2*n))/(trapets(2*n)-trapets(4*n));
    ratio_trap = [ratio_trap; ratio];
end

ratio_simp = [];
for n = series_of_n
    ratio = (simpson(n)-simpson(n*2))/(simpson(n*2)-simpson(n*4));
    ratio_simp = [ratio_simp; ratio];
end
disp(["Trapets" "Simpson"; 
			ratio_trap ratio_simp]);

%% 4c Trapetsregeln i 2D

trap_20 = trapets2d(1000);
disp("Volym = " + trap_20)

series_of_n = 60:2:80
ratio_trap = [];
for n = series_of_n 
    ratio = (trapets2d(n)-trapets2d(2*n))/(trapets2d(2*n)-trapets2d(4*n));
    ratio_trap = [ratio_trap; ratio];
end
disp(["Konvergensstudie"; ratio_trap])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V=trapets(n)
%  Indata:
%
%  n  - antal delintervall (skalär)
%
%  Utdata:
%
%  V - volymen (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	R = 3;
	g = @(r) 3* r^3 * exp(-r) / (1 + 1/3 * sin(8*r/5));

	spaced = linspace(0,R,n+1);
	h = R/(n);

	I = 0;
	for r = spaced(1:end-1)
				trap = h*(g(r+h)*(r+h)+g(r)*(r))/2;
				I = I + trap;
	end
	V_0 = g(R)*R^2*pi;
	V =  V_0 - 2*pi*I;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V=simpson(n)
%  Indata:
%
%  n  - antal delintervall (skalär)
%
%  Utdata:
%
%  V - volymen (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	R = 3;
	g = @(r) 3* r^3 * exp(-r) / (1 + 1/3 * sin(8*r/5));

	intervals = linspace(0,R,n+1);
	h = R/(n);

	I = g(0) * 0 + g(R)*R; 
	for j = 2:(size(intervals,2)-1)
		coeff = 2;
		if mod(j,2) == 0 
			coeff = 4;
		end
		
		I = I + coeff * g(intervals(j)) * intervals(j);
	end
	I = h/3 * I;
    
	V_0 = g(R)*R^2*pi;
	V =  V_0 - 2*pi*I;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function V=trapets2d(n)

%  Indata:
%
%  n  - antal delintervall i varje koordinatriktning (skalär)
%
%  Utdata:
%
%  V - volymen (skalär)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	R = 3;
	g = @(r) 3* r^3 * exp(-r) / (1 + 1/3 * sin(8*r/5));
	f = @(x,y) g(R) - g(sqrt(x^2+y^2));

	% Er kod här...
	L = 3 * sqrt(2);
	x_intervals = linspace(-L/2, L/2, n+1);
	y_intervals = linspace(-L/2, L/2, n+1);
	h = L/n;
	I = 0;
	for i = 1:length(x_intervals)
			for j = 1:length(y_intervals) 
					coeff = 1;
					if i == 1 || i == length(x_intervals) || j == 1 || j == length(y_intervals)
							coeff = 1/2;
							if (i == 1 && j == 1) || (i == 1 && j == length(y_intervals)) || (i == length(x_intervals) && j == 1) || (i == length(x_intervals) && j == length(y_intervals))
									coeff = 1/4;
							end
					end
					I = I + f(x_intervals(i),y_intervals(j))*coeff*h^2;
			end
	end

	V = I;
end


