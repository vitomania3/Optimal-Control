% NON LINEAR...
clear;
clc;

t0 = 0;
T = 5;
k = 1;
numOfPoints = 10;
xMin = 0;
xMax = 2;
L = 5;
alpha = 1;
epsilon = 5;
tol = 0.00001;

% Check x0 in [xMin, xMax].
if (alpha < xMin || alpha > xMax)
   error('x0 is out of range.'); 
end

% Check L > xMax.
if (L <= xMax)
   error('L must be greater than xMax.'); 
end

t = linspace(t0, T, numOfPoints);

% Optimal index
opt_index = 1;

% Array of functionals
u1_index = 0;
array_u1 = []; 

% Cell of graphics
    % Graphics structure
    % graphcs(1, :) - vector u1
    % graphcs(2, :) - vector u2
    % graphcs(3, :) - vector psi1
    % graphcs(4, :) - vector psi2
    % graphcs(5, :) - vector x1
    % graphcs(6, :) - vector x2
gph_index = 0;
graphics = {};

%% case, when u_1 >= 0, u_2 in [0, k], k > 0

% We will iterate on the segment
    % split_psi1 - splitting of the psi1.
f = @(T) (1 - exp(-T)) * (cosh(T) - 1) + sinh(T) * (1 - T - exp(-T));
if f(T) > 0
    split_psi1 = linspace(-(L - alpha) * sinh(T) / f(T), (epsilon * (cosh(T) - 1) - (L - alpha) * sinh(T)) / f(T), numOfPoints);
else
    split_psi1 = linspace((epsilon * (cosh(T) - 1) - (L - alpha) * sinh(T)) / f(T), -(L - alpha) * sinh(T) / f(T), numOfPoints);
end

for i = 1 : numOfPoints
    % psi2(0) = psi2_0 function.
    psi2_0 = (L - alpha + split_psi1(i) * (sinh(T) - T)) / (cosh(T) - 1);
   
    % Case, when psi2_0 > 0.
    if psi2_0 > 0
        % x1(t), x2(t), psi1(t), psi2(t) functions.
        x1 = @(t) split_psi1(i) * (t + exp(-t)) + (psi2_0 - split_psi1(i)) * cosh(t) + alpha - psi2_0;
        x2 = @(t) split_psi1(i) * (1 - exp(-t)) + (psi2_0 - split_psi1(i)) * sinh(t);
        psi2 = @(t) split_psi1(i) + (psi2_0 - split_psi1(i)) * exp(t);
        psi1 = @(t) split_psi1(i) * ones(1, numOfPoints);
        
        % u1(t), u2(t) functions.
        u1 = @(t) psi2(t);
        u2 = @(t) zeros(1, numOfPoints);
        
        % Case, when psi2_0 - psi1 > 0.
        if psi2_0 - split_psi1(i) > 0
            u1_index = u1_index + 1;
            array_u1(u1_index) = integral(@(t) u1(t) .* u1(t), t0, T);

            gph_index = gph_index + 1;
            graphics{gph_index} = [u1(t); u2(t); psi1(t); psi2(t); x1(t); x2(t)];
        end
    
        % Case, when psi2_0 = psi1.
        if abs(psi2_0 - split_psi1(i)) < tol
            if 0 <= (L - a) * (1 - exp(-T)) / (T + exp(-T) - 1) && ...
                    (L - a) * (1 - exp(-T)) / (T + exp(-T) - 1) <= epsilon
                
                psi2_0 = (L - alpha) / (T + exp(-T) - 1);
                psi2 = @(t) psi2_0 * ones(numOfPoiunts);
                x1 = @(t) psi2_0 * (1 - exp(-t));
                x2 = @(t) psi2_0 * (t + exp(-t) - 1) + alpha;

                u1_index = u1_index + 1;
                array_u1(u1_index) = integral(@(t) u1(t) .* u1(t), t0, T);

                gph_index = gph_index + 1;
                graphics{gph_index} = [u1(t); u2(t); psi1(t); psi2(t); x1(t); x2(t)];
            end
        end
        
        % Case, when psi2_0 - psi1 < 0.
        if (psi2_0 - split_psi1(i) < 0)
            % t_sw - time to switch
            t_sw = log(-split_psi1(i) / (psi2_0 - split_psi1(i)));
            
            new_psi2 = @(t) (split_psi1(i) / (1 + k)) * (1 - exp((1 + k) * (t - t_sw)));
            new_x2 = @(t) (split_psi1(i) * (1 - exp(-t_sw)) + (psi2_0 - split_psi1(i)) * sinh(t_sw)) * exp(-(1 + k) * (t - t_sw));
            F = @(t) (1 / (1 + k)) * (exp(t_sw) - 1 - sinh(t_sw)) * (1 - exp(-(1 + k) * (t - t_sw))) - ...
                cosh(t_sw) + exp(t_sw) * (t_sw - 1) + 2;
            new_x1 = @(t) (psi2_0 / (exp(t_sw) - 1)) * F(t) + alpha;
            psi2_0 = (L - alpha) * (exp(t_sw) - 1) / F(T);

            vecOfZeros = zeros(1, numOfPoints - size(t(t < t_sw), 2));
            vecOfOnes = ones(1, numOfPoints - size(t(t < t_sw), 2));
            
            u1_index = u1_index + 1;
            array_u1(u1_index) = trapz(t, cat(2, u1(t(t < t_sw)), vecOfZeros).^2);
            
            gph_index = gph_index + 1;
            graphics{gph_index} = [cat(2, u1(t(t < t_sw)), vecOfZeros);
                                   cat(2, zeros(1, size(t(t < t_sw), 2)), k * vecOfOnes);
                                   psi1(t);
                                   cat(2, psi2(t(t < t_sw)), new_psi2(t(t >= t_sw))); 
                                   cat(2, x1(t(t < t_sw)), new_x1(t(t >= t_sw)));
                                   cat(2, x2(t(t < t_sw)), new_x2(t(t >= t_sw)))];
        end
    end
end

% Case, when psi2_0 = 0.
if 0 <= (L - alpha) * (1 - cosh(T)) / (T - sinh(T)) && ...
        (L - alpha) * (1 - cosh(T)) / (T - sinh(T)) <= epsilon
    psi1 = (L - alpha) / (T - sinh(T));
    psi2 = @(t) psi1 * (1 - exp(t));
    x2 = @(t) psi1 * (1 - cosh(t));
    x1 = @(t) psi1 * (t - sinh(t)) + alpha;

    u1_index = u1_index + 1;
    u1 = @(t) psi2(t);
    u2 = @(t) zeros(1, numOfPoints);
    array_u1(u1_index) = integral(@(t) u1(t) .* u1(t), t0, T);

    gph_index = gph_index + 1;
    graphics{gph_index} = [u1(t); u2(t); psi1 * ones(1, numOfPoints); psi2(t); x1(t); x2(t)];
end
    
% Case, when psi2_0 < 0.
suppFunc = @(t_st) x2(T, t_st);
tmpSplit = linspace(t0, T, 2 * numOfPoints);
for j = 1 : 2 * numOfPoints - 1
    if suppFunc(tmpSplit(j)) < epsilon
        t_st = tmpSplit(j);
        psi1 = ((L - alpha) * (1 + k)^3) / ((T - t_st) * (1 + k) - sinh((1 + k) * (T - t_st)));
        psi2 = @(t) (psi1(t_st) / (1 + k)) * (1 - exp((1 + k) * (t - t_st)));
        x2 = @(t) (psi1(t_st) / (1 + k)^2) * (1 - cosh((1 + k) * (t - t_st)));
        x1 = @(t) (psi1(t_st) / (1 + k)^2) * (t - t_st - (1 / (1 + k)) * sinh((1 + k) * (t - t_st))) + alpha; 
        
        u1_index = u1_index + 1;
        u1 = @(t) psi2(t);
        u2 = @(t) zeros(1, numOfPoints);
        array_u1(u1_index) = integral(@(t) u1(t(t >= t_st)) .* u1(t(t >= t_st)), t0, T);

        gph_index = gph_index + 1;
        graphics{gph_index} = [cat(2, zeros(1, size(t(t < t_st), 2)), u1(t(t >= t_st))); 
                               u2((t >= t_st)); 
                               psi1 * ones(1, 
                               size((t >= t_st), 2)); 
                               psi2(t(t >= t_st)); 
                               x1(t(t >= t_st)); 
                               x2(t(t >= t_st))];
        break;
    end
end

% Search of the minimum value of the functional
tmp = array_u1(1);
for i = 2 : size(array_u1, 2)
    if array_u1(i) < tmp
        tmp = array_u1(i);
        opt_index = i;
    end
end

% Derivation of the minimum value of the functional
array_u1(opt_index)

plot(graphics{opt_index}(6, :), graphics{opt_index}(4, :));
grid on;
%%

% We will iterate on the segment
f = @(T) (1 - exp(-T)) * (cosh(T) - 1) + sinh(T) * (1 - T - exp(-T));
if f(T) > 0
    split_psi1 = linspace(-(L - alpha) * sinh(T) / f(T), (epsilon * (cosh(T) - 1) - (L - alpha) * sinh(T)) / f(T), numOfPoints);
else
    split_psi1 = linspace((epsilon * (cosh(T) - 1) - (L - alpha) * sinh(T)) / f(T), -(L - alpha) * sinh(T) / f(T), numOfPoints);
end

psi2_0 = (L - alpha + split_psi1 * (sinh(T) - T)) / (cosh(T) - 1);
x1 = @(t, i) split_psi1(i) * (t + exp(-t)) + (psi2_0(i) - split_psi1(i)) * cosh(t) + alpha - psi2_0(i);
x2 = @(t, i) split_psi1(i) * (1 - exp(-t)) + (psi2_0(i) - split_psi1(i)) * sinh(t);
psi2 = @(t, i) split_psi1(i) + (psi2_0(i) - split_psi1(i)) * exp(t);

for i = 1 : numOfPoints
    % case, when psi2_0 = split_psi1.
    if (abs(psi2_0(i) - split_psi1(i)) < tol) && (0 <= ((L - alpha) * (1 - exp(-T))) / (T + exp(-T) - 1) && ...
        ((L - alpha) * (1 - exp(-T))) / (T + exp(-T) - 1) <= epsilon)
        psi2_0 = (L - alpha) / (T + exp(-T) - 1);
        psi2 = @(t) psi2_0 * ones(1, size(t, 2));
        x1 = @(t) psi2_0 * (1 - exp(-t));
        x2 = @(t) psi2_0 * (t + exp(-t) - 1) + alpha;
    %     plot(x1(t), x2(t));
        plot(x2(t), psi2(t));
    elseif (psi2_0(i) - split_psi1(i) < 0)
        t_sw = log(-split_psi1(i) / (psi2_0(i) - split_psi1(i)));

        ps1 = split_psi1(i);
        odeFunc = @(t,y) [y(2) - y(1); ps1+ y(2)];
        [tt, xx] = ode45(odeFunc, linspace(t0, T, numOfPoints), [ps1 * (1 - exp(-t_sw)) + (psi2_0(i) - ps1) * sinh(t_sw); 0]);

        plot(xx(:, 1), xx(:, 2), 'r');
        hold on;
        plot(x2(t_sw, i), psi2(t_sw, i), 'Og');
    elseif (psi2_0(i) - split_psi1(i) > 0)
        plot(x2(t, i), psi2(t, i), 'b');
    elseif (abs(psi2_0) < tol) && (0 <= (L - alpha) * (1 - cosh(T)) / (T - sinh(T)) && ...
           (L - alpha) * (1 - cosh(T)) / (T - sinh(T)) <= epsilon)
        split_psi1 = (L - alpha) / (T - sinh(T))
        psi2 = @(t) split_psi1 * (1 - exp(t));
        x2 = @(t) split_psi1 * (1 - cosh(t));
        x1 = @(t) split_psi1 * (t - sinh(t)) + alpha;
        plot(x2(t), psi2(t), 'g');
    %     plot(x1(t), x2(t));
    end
    hold on;
end
grid on;

%% 
epsilon = 10;
split_psi1 = @(t_st) ((L - alpha) * (1 + k)^3) ./ ((T - t_st) * (1 + k) - sinh((1 + k) .* (T - t_st)));
x2 = @(t_st) (split_psi1(t_st) ./ (1 + k)^2) .* (1 - cosh((1 + k) .* (T - t_st))) - epsilon;
fzero(x2, [t(2), t(size(t, 2) - 1)]);

%% case, when u_1 in R, u_2 in [0, k], k > 0