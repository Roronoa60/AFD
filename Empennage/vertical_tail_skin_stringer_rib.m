clear
clc
addpath('interpolation_functions')

%multi_cell_shear_flow(0, 'VT');

%[weights, best_geometry] = main_optimisation(78, 'HT')

weights_array = D_section_optimisation(78, 265, 'HT', 'top');

function [weights_array] = D_section_optimisation(E, tau_tres, VT_or_HT, load_direction)
    % Assume semicircular D_section
    wing_span = extract_dimension(0, 'wing_span', VT_or_HT);
    weights_array = {};
    
    for num_rib = 7 % 0:20
        num_rib
        cell_span = wing_span/(num_rib + 1); % a
        t_dist = [];
        total_weight = 0;
        
        for section = 0:(num_rib)
            span_location = section*cell_span;
            
            LE_span = extract_dimension(span_location, 'chord', VT_or_HT)*0.1; % 10% LE Assumption
            R = LE_span*1; % D cell radius of curvature (guess)
            
            bh = extract_dimension(span_location, 'bh', VT_or_HT);
            
            t = 5; % D cell Initial Estimate thickness (mm)
            
            % pi/2 / rt(2) ~= 1.11
            b = 1.1*sqrt((bh/2)^2 + LE_span^2); % curve_length
            
            for i = 1:6
                b_Rt = b/sqrt(R*t/1000);
                a_b = cell_span/b;

                Ks = extr_nose_k(b_Rt,a_b); % a=a/b and b= b/sqrt(Rt)
                t = sqrt((tau_tres*10^6*b^2)/(Ks*E*10^9))*1000;
            end
            if t < 1
                t = 1;
            end
            t_dist(end+1) = t;
            rib_weight = (2/3)*bh*LE_span*1/1000; % 2/3 numerical estimate
            skin_weight = cell_span*b*2*t/1000;
            total_weight = total_weight + rib_weight + skin_weight;
        end
            
        weights_array(end + 1, :) = {num_rib, cell_span, total_weight, t_dist};
        
                
    end
end

function [t,sigma] = simple_panel_buckling(M, c, bh, E)
    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    t = (N/(3.62*E*10^9)*c^2)^(1/3)*1000; % t - Skin Thickness (mm)
    sigma = N/t/1000; % sigma - Buckling stress (N/mm^2)
end

function weight = stringer_panel_weight(As, c, n, t)
    weight = n*As + c*t;
end

function [weights, best_geometry] = main_optimisation(E, VT_or_HT)
    wing_span = extract_dimension(0, 'wing_span', VT_or_HT);
    bh_min = extract_dimension(wing_span, 'bh', VT_or_HT); %bh_min coresponds to box height at wing tip
    
    % Find shear flow, box height, box width and bending moment at root (i.e. 0 span)
    q = skin_shear_flow(0, VT_or_HT);
    bh = extract_dimension(0, 'bh', VT_or_HT);
    c = extract_dimension(0, 'c', VT_or_HT);
    M = extract_force(0, 'BM', VT_or_HT, 'top');
    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    
    % Simple Panel Weight
    [t_simple, sigma_simple] = simple_panel_buckling(M, c, bh, E); % returns skin thickness (mm) and buckling stress (N/mm^2) for simple panel
    simple_panel_weight = c*t_simple*1000;

    weights = {}; % {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'root weight', 'total weight'}
    titles = {'n', 't', 'Rt', 'Rdh', 'Rb', 'F', 'L', 'root weight', 'total weight', 'As/bt', 't_dist', 'L_dist', 'Tr_dist', 'te'};
    
    best_geometry = struct('total_weight', 1000);
        
    for n = 20
        n
        for Rt = [0.5,0.6,0.7,0.8,0.9,1,1.25,1.5,2]
            for Rdh = [0.3,0.4,0.5]
                for Rb = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
                    
                    b = c/(n+1); % b - stringer pitch (m)
                    h = Rb*b; % h - stringer height (m)
                    d = Rdh*h; % d - flange height (m)
                    
                    if (h < bh_min/2) & (2*d < b) % Check stringer height within min box height
                        % Get Optimum FARRAR
                        As_bt = Rt*(h+2*d)/b; % As/(b*t) = ((h + 2*d)*ts)/(b*t) = Rt*(h+2*d)/b
                        F = Farrar_Calc(As_bt, Rt);
                        if F > 0.6
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

                            [L_distribution, Tr_distribution, t_distribution] = skin_rib_distribution(geometry, E, VT_or_HT);
                            [rib_weight, stringer_weight, total_weight, L] = tail_weight(geometry, L_distribution, Tr_distribution, t_distribution, VT_or_HT);
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
    eqn = (q/t / (Ks * E * (t/g.b)^2))^2 + N/( t*(1 + Rt*Rb) ) / (Kc * E * (t/g.b)^2 ) == 0.99;
    S = solve(eqn, t, 'Real', true); % Solve for t
    root = 1000*eval(min(S(S>0))); % Find smallest real positive solution for t
    if isempty(root)
        t = 1000;
    else
        t = root;
    end        
end

function [L_distribution, Tr_distribution, t_distribution] = skin_rib_distribution(geometry, E, VT_or_HT) % L varying along span
    g = geometry;
    F = geometry.F;
    
    wing_span = extract_dimension(0, 'wing_span', VT_or_HT);
        
    Tr_array = [];
    t_array = [g.t];
    L_distribution = [0];
    
    bh = extract_dimension(0, 'bh', VT_or_HT);
    c = extract_dimension(0, 'c', VT_or_HT);
    M = extract_force(0, 'BM', VT_or_HT, 'top');
    N = M/(c*bh);
    while 1 % iterate through each section in wing and assign rib spacing, rib thickness and skin thickness
        Kc = extr_zstr_k(g.h/g.b, g.ts/t_array(end), g.d/g.h);
        te = t_array(end) + g.As/(1000*g.b);
        L = 10^3*(3.62/Kc)^2*F^2*te^2*E/N; % Solution to Farrar Equation
        if abs(L_distribution(end) + L) >= wing_span % check if section beyond wing_span
            break
        else
            span_location = L_distribution(end) + L;
            q = skin_shear_flow(span_location, VT_or_HT);
            M = extract_force(span_location, 'BM', VT_or_HT, 'top');
            c = extract_dimension(span_location, 'c', VT_or_HT);
            bh = extract_dimension(span_location, 'bh', VT_or_HT);
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

function [rib_weight, stringer_weight, total_weight, L] = tail_weight(g, L_distribution, Tr_distribution, t_distribution, VT_or_HT)
    rib_weight = 0;
    stringer_weight = 0;
    total_weight = 0;
    wing_span = extract_dimension(0, 'wing_span', VT_or_HT);
    L_distribution(end + 1) = wing_span;
    
    for i = 1:(size(L_distribution,2)-1)
        span_start = L_distribution(i);
        span_end = L_distribution(i+1);
        section_length = L_distribution(i+1)-L_distribution(i);
        % bh1 = extract_dimension(span_start, 'bh', VT_or_HT);
        bh = extract_dimension(span_end, 'bh', VT_or_HT);
        c1 = extract_dimension(span_start, 'c', VT_or_HT);
        c2 = extract_dimension(span_end, 'c', VT_or_HT);
        
        te = t_distribution(i) + g.As/(1000*g.b);
        
        stringer_section_weight = (te/1000 *section_length*(c1+c2)/2); 
        if VT_or_HT == 'VT'
            stringer_section_weight = stringer_section_weight*2; % Multiply by 2 for vertical tail
        end
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

function q = skin_shear_flow(L_span, VT_or_HT) % returns q (N/mm)
% shear flow is purely from torque on the top and bottom skin
    T = extract_force(L_span, 'T', VT_or_HT, 'top');
	c = extract_dimension(L_span, 'c', VT_or_HT);
    b2 = extract_dimension(L_span, 'bh', VT_or_HT);
    q = T/(2*b2*c*1000);  
end

%{
function q = multi_cell_shear_flow(L_span, VT_or_HT)
    syms t real
    
    T = extract_force(L_span, 'T', VT_or_HT, 'top')*10^-3;
    SF = extract_force(L_span, 'SF', VT_or_HT, 'top');
    c = extract_dimension(L_span, 'c', VT_or_HT);
    
    hw = extract_dimension(L_span, 'bh', VT_or_HT);
    LE_span = extract_dimension(L_span, 'chord', VT_or_HT)*0.1; % 10% LE Assumption
    
    Sn = 1.1*sqrt((hw/2)^2 + LE_span^2)*2*10^3;
    Sr = 2*c*10^3; % top and bottom lengths combined
    
    An = (2/3)*hw*LE_span*10^6;
    Ar = hw*c*10^6;
    
    qw_front = 80; % Place holder
    qw_rear = 40; %Place holder
    tw_front = 3; %Place holder
    tw_rear = 3; %Place holder
    
    tn = 3; %Place holder
    
    matrix = [2*An      , Ar     ;
              Sn/(An*tn), -Sr/(Ar*t) ];
    RHS = [T-Ar*qw_rear ;
           -(1/Ar)*(hw*qw_rear/tw_rear)-(1/Ar + 1/An)*(hw*qw_front/tw_front)]
    solution = matrix'*RHS
    qn = eval(solution(1))
    qr(t) = eval(solution(2));
    eval(qr(3))
    
end
%}

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