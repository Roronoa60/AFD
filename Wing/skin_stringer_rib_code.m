
clear
clc

[weights_array, best_geometry] = main_optimisation(78, 30);
addpath('Stringer Panel Result Data')
load('n30_best_geometry.mat')
[weights, best_n] = bottom_stringer_panel(78, 30, best_geometry, best_geometry.L_distribution);


function  [t,sigma] = simple_panel_buckling(M, c, bh, E)
    % M - Bending Moment at section (Nm)
    % c - width of box (m)
    % bh - Web (box) height (m) 
    % E - Young Modulus of panel (Pa)
    
    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    % Assuming that bending about bottom surface to be conservative so moment arm = b2 
    
    t = (N/(3.62*E*10^9)*c^2)^(1/3)*1000; % t - Skin Thickness (mm)
    sigma = N/t/1000; % sigma - Buckling stress (N/mm^2)
end

function  [weights, best_geometry] = main_optimisation(E , num_stringers)
    % M - root bending moment (Nm)
    % c - root box width (m)
    % bh - root box height (m)
    % bh_min - min box height
    % E - Young Modulus of panel and stringer material (Pa)
    % q - Shear Flow

    wing_span = extract_dimension(0, 'wing_span');
    bh_min = extract_dimension(wing_span, 'bh'); %bh_min coresponds to box height at wing tip
    
    % Find shear flow, box height, box width and bending moment at root (i.e. 0 span)
    q = skin_shear_flow(0);
    bh = extract_dimension(0, 'bh');
    c = extract_dimension(0, 'c');
    M = extract_force(0, 'BM', 'top');
    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    
    % Simple Panel Weight
    [t_simple, sigma_simple] = simple_panel_buckling(M, c, bh, E); % returns skin thickness (mm) and buckling stress (N/mm^2) for simple panel
    simple_panel_weight = c*t_simple*1000;
       
    % n - No. of stringers
    % Rt = ts/t = stringer web thickness/skin thickness
    % Rdh = d/h = flange width/flange height
    % Rb = d*/b = h/b stringer height/stringer pitch  (d* as Gp 17 denote stringer height as d not h)
    
    weights = {}; % {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'root weight', 'total weight'}
    titles = {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'L', 'root weight', 'total weight', 'As/bt', 't_dist', 'L_dist', 'Tr_dist', 'te'};
    
    best_geometry = struct('total_weight', 1000);
        
    for n = num_stringers %30 %6:2:50
        n
        for Rt = [0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,2]
            for Rdh = [0.3,0.4,0.5]
                for Rb = 0.5 % [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5]
                    
                    b = c/(n+1); % b - stringer pitch (m)
                    h = Rb*b; % h - stringer height (m)
                    d = Rdh*h; % d - flange height (m)
                    
                    if (h < bh_min/2) & (2*d < b) % Check stringer height within min box height
                        % Get Optimum FARRAR
                        As_bt = Rt*(h+2*d)/b; % As/(b*t) = ((h + 2*d)*ts)/(b*t) = Rt*(h+2*d)/b
                        F = Farrar_Calc(As_bt, Rt);
                        if F > 0.5
                            geometry = struct('n', n, 'b', b, 'Rt', Rt, 'h', h, 'd', d);
                            % Solve for skin thickness using combined loads inequality
                            t = solve_skin_thickness(geometry, c, E, 0.5, N, q, 'root');
                            ts = Rt*t;
                            
                            if ts < 1
                                ts = 1;
                                t = ts/Rt;
                            end

                            % Calculate Weight Saving
                            As = As_bt * t * b * 1000;
                            stringer_panel_weight = n*As + c*t*1000; % Calculates Stringer Panel Weight
                            root_weight_save = (simple_panel_weight - stringer_panel_weight)/simple_panel_weight;

                            te = t + As/(1000*b);
                            geometry = struct('n', n, 'b', b, 't', t, 'ts', ts, 'te', te, 'h', h, 'd', d, 'F', F, 'As', As, 'weight', root_weight_save, 'Rt', Rt, 'Rdh', Rdh, 'Rb', Rb);

                            [L_distribution, Tr_distribution, t_distribution] = varying_L_solution(geometry, E);
                            [rib_weight, stringer_weight, total_weight, L] = wing_weight(geometry, L_distribution, Tr_distribution, t_distribution);
                            weights(end + 1, :) = {n, t, Rt, Rdh, Rb, F, L, root_weight_save, total_weight*1000, As/(b*t*1000), t_distribution, L_distribution, Tr_distribution, te}; 
                            if total_weight*1000 < best_geometry.total_weight
                                best_geometry = geometry;
                                best_geometry.total_weight = total_weight*1000;
                                best_geometry.As_bt = As/(b*t*1000);
                                best_geometry.t_distribution = t_distribution;
                                best_geometry.L_distribution = L_distribution;
                                best_geometry.Tr_distribution = Tr_distribution;
                            end
                        end
                    end
                end
            end
            Rt
            size(weights,1)
        end
    end
    
end

function [weights, best_n] = bottom_stringer_panel(E, num_stringers, best_geometry, L_distribution)
    bg = best_geometry;
    weights = {};
    best_n = struct('n', 0, 'stringer_weight', 10000); % dummy default
    
    c = extract_dimension(0, 'c');
    for n = 6:2:40
        n
        c = extract_dimension(0, 'c');
        bg.b = c/(n+1);
        bh = extract_dimension(0, 'bh');
        M = extract_force(0, 'BM', 'bottom');
        N = M/(c*bh);
        q = skin_shear_flow(0);
        
        tn = bg.t;
        for i = 1:8
            tn = solve_skin_thickness(bg, c, E, L_distribution(2), abs(N), q, 'span', tn);
        end
        if tn < 1            tn = 1;
        end
        t_array = [tn];

        for section = 2:size(L_distribution,2)
            L = L_distribution(section) - L_distribution(section - 1);
            span_location = L_distribution(section);
            q = skin_shear_flow(span_location);
            M = extract_force(span_location, 'BM','bottom');
            c = extract_dimension(span_location, 'c');
            bh = extract_dimension(span_location, 'bh');
            N = M/(c*bh); % compressive Load at section*L
            t_array(end + 1) = solve_skin_thickness(bg, c, E, L, abs(N), q, 'span', t_array(end));
            if t_array(end) < 1
                t_array(end) = 1;
            end
        end
        [rib_weight, stringer_weight, total_weight, L] = wing_weight(bg, L_distribution, t_array, t_array);
        weights(end + 1, :) = {n, t_array(1), bg.F, stringer_weight*1000, t_array};
        if stringer_weight*1000 < best_n.stringer_weight
            best_n.n = n;
            best_n.stringer_weight = stringer_weight*1000;
            best_n.t_distribution = t_array;
        end
    end
end

function t = solve_skin_thickness(g, c, E, L, N, q, type, t_skin) % Use combined loading inequality of compressive and shear stress to solve for t % NOT TESTED 
    if type == 'root'
        Rt = g.Rt;
    else
        Rt = g.ts/t_skin;
    end
        
    E = E*10^9;
   
    syms t
    Rb = g.h / g.b;
    Kc = extr_zstr_k(g.h/g.b, Rt, g.d/g.h);
    Ks = extr_ks( L/c, 1);
    eqn = (q/t / (Ks * E * (t/g.b)^2))^2 + N/( t*(1 + Rt*Rb) ) / (Kc * E * (t/g.b)^2 ) == 0.99;    S = solve(eqn, t, 'Real', true); % Solve for t
    root = 1000*eval(min(S(S>0))); % Find smallest real positive solution for t
    if isempty(root)
        t = 1000;
    else
        t = root;
    end        
end

function [L_distribution, Tr_distribution, t_distribution] = varying_L_solution(geometry, E) % L varying along span
    g = geometry;
    F = geometry.F;
    
    wing_span = extract_dimension(0, 'wing_span');
        
    Tr_array = [];
    t_array = [g.t];
    L_distribution = [0];
    
    bh = extract_dimension(0, 'bh');
    c = extract_dimension(0, 'c');
    M = extract_force(0, 'BM', 'top');
    N = M/(c*bh);
    while 1 % iterate through each section in wing and assign rib spacing, rib thickness and skin thickness
        Kc = extr_zstr_k(g.h/g.b, g.ts/t_array(end), g.d/g.h);
        te = t_array(end) + g.As/(1000*g.b);
        L = 10^3*(3.62/Kc)^2*F^2*te^2*E/N; % Solution to Farrar Equation
        if L_distribution(end) + L > wing_span % check if section beyond wing_span
            break
        else
            span_location = L_distribution(end) + L;
            q = skin_shear_flow(span_location);
            M = extract_force(span_location, 'BM','top');
            c = extract_dimension(span_location, 'c');
            bh = extract_dimension(span_location, 'bh');
            N = M/(c*bh); % compressive Load at section*L
            
            L_distribution(end + 1) = span_location;
            Tr_array(end + 1) = rib_thickness(geometry, c, E, bh, L, M, t_array(end));
            t_array(end + 1) = solve_skin_thickness(geometry, c, E, L, N, q, 'span', t_array(end));
            if t_array(end) < 1
                t_array(end) = 1;
            end
            if Tr_array(end) < 1
                Tr_array(end) = 1;
            end
        end
    end

    Tr_distribution = Tr_array;
    t_distribution = t_array;
end

function [L_distribution, Tr_distribution, t_distribution] = constant_L_solution(geometry, E) % L constant along spa
    g = geometry;
    F = geometry.F;
    
    wing_span = extract_dimension(0, 'wing_span');
    L_distribution = [0];
    Tr_array = [];
    t_array = [g.t];
    
    
    bh = extract_dimension(0, 'bh');
    c = extract_dimension(0, 'c');
    M = extract_force(0, 'BM', 'top');
    N = M/(c*bh);
    D = bh;
    L = 0.5;
    
    for i = 1:20
        Kc = extr_zstr_k(g.h/g.b, g.ts/t_array(end), g.d/g.h);
        Tr = rib_thickness(geometry, c, E, bh, L, M, t_array(end))/1000;
        
        L = (4*F^2*D^2*Tr^2*E*10^9/N)^(1/3);
        te = 1000*sqrt((N*(Kc/3.62)^2*L)/(F^2*E*10^9));
        t = te - g.As/(1000*g.b);
        t_array = [t];
    end
    
    if t_array(end) < g.t
        t_array(end) = g.t;
    end
        
    while L_distribution(end) + L < wing_span % iterate through each section in wing and assign rib spacing, rib thickness and skin thickness
        span_location = L_distribution(end) + L;
        q = skin_shear_flow(span_location);
        M = extract_force(span_location, 'BM', 'top');
        c = extract_dimension(span_location, 'c');
        bh = extract_dimension(span_location, 'bh');
        N = M/(c*bh); % compressive Load at section*L
        Kc = extr_zstr_k(g.h/g.b, g.ts/t_array(end), g.d/g.h);
        
        L_distribution(end + 1) = span_location;
        Tr_array(end + 1) = rib_thickness(geometry, c, E, bh, L, M, t_array(end));
        te = 1000*sqrt((N*(Kc/3.62)^2*L)/(F^2*E*10^9));
        t1 = te - g.As/(1000*g.b);
        t2 = solve_skin_thickness(geometry, c, E, L, N, q, 'span', t_array(end));
        if t1 > t2
            t_array(end + 1) = t1; 
        else
            t_array(end + 1) = t2;
        end
        if t_array(end) < 1
            t_array(end) = 1;
        end
        if Tr_array(end) < 1
            Tr_array(end) = 1;
        end
    end
    t_array(end) = 0;
    Tr_distribution = Tr_array;
    t_distribution = t_array;
end

function [rib_weight, stringer_weight, total_weight, L] = wing_weight(g, L_distribution, Tr_distribution, t_distribution)
    rib_weight = 0;
    stringer_weight = 0;
    total_weight = 0;
    wing_span = extract_dimension(0, 'wing_span');
    L_distribution(end + 1) = wing_span;
    
    for i = 1:(size(L_distribution,2)-1)
        span_start = L_distribution(i);
        span_end = L_distribution(i+1);
        section_length = L_distribution(i+1)-L_distribution(i);
        bh = extract_dimension(span_end, 'bh');
        c1 = extract_dimension(span_start, 'c');
        c2 = extract_dimension(span_end, 'c');
        
        te = t_distribution(i) + g.As/(1000*g.b);
        
        stringer_section_weight = (te/1000 *section_length*(c1+c2)/2); 

        if i == size(L_distribution,2)-1
            rib_section_weight = 0;
        else
            rib_section_weight = Tr_distribution(i)*c2*bh/1000;
        end
        stringer_weight = stringer_weight + stringer_section_weight;
        rib_weight = rib_weight + rib_section_weight;
        total_weight = total_weight + stringer_section_weight + rib_section_weight;      
    end
    L = L_distribution(2);
end

function q = skin_shear_flow(L_span) % returns q (N/mm)
% shear flow is purely from torque on the top and bottom skin
    T = extract_force(L_span,'T', 'top');
	c = extract_dimension(L_span, 'c');
    b2 = extract_dimension(L_span, 'bh');
    q = T/(2*b2*c*1000);  
end

function I = second_moment_area(c, te, hc) % c:(m),te(mm),hc(mm),I(m^4)
    % hc = beam depth between panels = box height
    I = c*(te/1000)^3/12 + c*(te/1000)*(hc/2)^2;
end

function Tr = rib_thickness(g, c, E, hc, L, M, skin_thickness) % Returns Tr in (mm)
% Calculates the rib thickness based on parameters
    % E in terms of GPa
    % hc is height of the wing box (hb)
    s = L;
    E = E*10^9;

    te = skin_thickness + g.As/(1000*g.b);
    I = second_moment_area(c, te, hc);  % Second moment of area of the skin
    F_crush = (M^2*s*hc*te/1000*c)/(2*E*I^2);
    Tr = ((F_crush*hc^2)/(3.62*c*E))^(1/3)*1000;
end