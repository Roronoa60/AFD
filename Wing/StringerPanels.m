clear
clc
close all

[t, sigma] = simple_panel_buckling(23286.8, 1, 0.2, 70E9);

%stringer_panel_buckling_thickness(10, 1, 1, 50, 23286.8, 1, 0.2, 70E9);

%{
Rib and Second moment of area Function Tests:
geometry = struct('t', 2, 'As', 0, 'b', 1);
c = 1;
E = 72;
hc = 0.2;
L = 0.5;
M = 23286.8;
rib_thickness(geometry, c, E, hc, L, M) % should equal 0.661 mm
%}

geometry = struct('n', 10, 'b', 345, 't', 15, 'ts', 10, 'h', 50, 'd', 15, 'F', 0.82, 'As', 80, 'weight', 0.611);
E = 70;

[L, Tr_distribution, t_distribution] = varying_L_solution(geometry, E);

scatter(L(2:end), Tr_distribution);
ylabel('Rib Thickness (mm)')
scatter(L(1:end), t_distribution);
xlabel('Wing Span (m)')
ylabel('skin Thickness (mm)')
ylim([0 16])

%% 3.1
function  [t,sigma] = simple_panel_buckling(M, c, bh, E)
    % M - Bending Moment at section (Nm)
    % c - width of box (m)
    % bh - Web (box) height (m) 
    % E - Young Modulus of panel (Pa)
    
    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    % Assuming that bending about bottom surface to be conservative so moment arm = b2 
    
    t = (N/(3.62*E)*c^2)^(1/3)*1000; % t - Skin Thickness (mm)
    sigma = N/t/1000; % sigma - Buckling stress (N/mm^2)
end

%% 3.2.1 - 3.2.5
function  [min_weight_geometry, max_F_geometry] = root_geometry(M, c, bh, bh_min, E, q) %(n, t2, ts, h, M, c, b2, E)
    % M - root bending moment (Nm)
    % c - root box width (m)
    % bh - root box height (m)
    % bh_min - min box height
    % E - Young Modulus of panel and stringer material (Pa)
    % q - Shear Flow

    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    
    %% Simple Panel Weight
    [t_simple, sigma_simple] = simple_panel_buckling(M, c, bh, E); % returns skin thickness and buckling stress for simple panel
    simple_panel_weight = c*t_simple*1000;
    
    %%    
    % weights = {'n', 'b', 't', 'ts', 'h', 'd', 'F', 'As', 'weight'};
    min_w = struct('n', 0, 'b', 0, 't', 0, 'ts', 0, 'h', 0, 'd', 0, 'F', 0, 'As', 0, 'weight', 0); % minimum weight geometry
    max_F = struct('n', 0, 'b', 0, 't', 0, 'ts', 0, 'h', 0, 'd', 0, 'F', 0, 'As', 0, 'weight', 0); % maximum F geometry
    
    % n - No. of stringers
    % Rt = ts/t = stringer web thickness/skin thickness
    % Rdh = d/h = flange width/flange height
   
    % Rb = d*/b = h/b stringer height/stringer pitch  (d* as Gp 17 denote stringer height as d not h)
    
    for n = 1:16
        for Rt = [0.5,0.6,0.7,0.9,1,1.25,1.5,2]
            for Rdh = [0.3,0.4,0.5]
                for Rb = [0.05, 0.1,0.15,0.2,0.25,0.3]
                    
                    b = c/(n+1); % b - stringer pitch
                    h = Rb*b; % h - stringer height
                    d = Rdh*h; % d - flange height
                    
                    if h < bh_min/2 % Check stringer height within min box height
                        
                        % Get Optimum FARRAR
                        As_bt = Rt*(h+2*d)/b; % As/(b*t) = ((h + 2*d)*ts)/(b*t) = Rt*(h+2*d)/b
                        F = FARRAR(Rt, As_bt); %!!!!!!!!!!!!!
                        
                        % Solve for skin thickness using combined loads inequality
                        t = solve_skin_thickness_1(N, Rt, Rb, b, E, q);
                        ts = Rt*t;
                        
                        % Calculate Weight Saving
                        As = As_bt * t * b;
                        stringer_panel_weight = stringer_panel_weight(n, b, t, As);
                        weight_save = (simple_panel_weight - stringer_panel_weight)/simple_panel_weight;
                        
                        % If new geometry improves weight saving
                        if min_w.weight < weight_save
                            min_w = struct('n', n, 'b', b, 't', t, 'ts', ts, 'h', h, 'd', d, 'F', F, 'As', As, 'weight', weight_save);
                        end
                        % If new geometry improves F
                        if max_F.F < F | (max_F.F == F & max_F.weight > weight_save) % if F matches, pick lighter geometry
                            max_F = struct('n', n, 'b', b, 't', t, 'ts', ts, 'h', h, 'd', d, 'F', F, 'As', As, 'weight', weight_save);
                        end
                        % weights(end+1,:) = {n, b, t, ts, h, d, F, weight_save}
                    else
                        % print out some message that h is too large?
                    end
                end
            end
        end
    end
    min_weight_geometry = min_w;
    max_F_geometry = max_F;
    % Purpose of storing all the variables in a cell array is the option to store all different iterations in one array
    %{
    n = min_weight(1);
    b = min_weight(2);
    t = min_weight(3);
    ts = min_weight(4);
    h = min_weight(5);
    d = min_weight(6);
    F = min_weight(7);
    As = min_weight(8);
    %}
    %{
    % t2 - Skin Thickness (mm)
    % ts - Stringer web thickness (mm)
    % h - Depth of stringer (mm)
    % M - Bending Moment at section (Nm)
    % c - Box width (m)
    % b2 - Web height (m)
    % E - Young Modulus of panel (Pa)
    
    % b - Stringer spacing (mm)
    b = c/n*1000;
    % d - Flange width of stringer (mm)
    d = 0.3*h;
    % As- Cross-section area of stringer (mm^2)
    As = ts*h + 2*d*ts;
    % t_e - Effective thickness (mm)
    t_e = t2 + As/b;
    % b1 - stringer pitch (m)
    b1 = c/(n+1);
    
    %[buckling_thickness, buckling_stress] = simple_panel_buckling(M, b1, b2, E)
    N = M/(c*b2); % N - Compressive Load of unit length (N/m)
    t_b = (N/(3.62*E)*b1^2)^(1/3)*1000; % t_b - Buckling Skin Thickness (mm)
    sigma_0 = N/t_b/1000; % sigma - Buckling stress (N/mm^2)
    
    ts/t2; % ts/t2
    As/(b*t2); % As/bt
    sigma_ratio = 1.38; % sigma_ratio = sigma_cr/sigma_0 = Catchpole_fun(ts/t2,As/(b*t2))   READ FROM CATCHPOLE DIAGRAM
     
    % FARRAR efficiency factor
    F = 0.82; % F = FARRAR_Factor(ts/t2, As/(b*t2))
    % sigma_f = F/sqrt(L/N/E)/10^6 % sigma_f - Mean stress by skin and stringer
    sigma_cr = sigma_0*sigma_ratio; % sigma_cr - critical buckling stress
    
    % Equate critical initial buckling stress with critical flexural
    % buckling and rearange for L:
    L = F^2*N*E/(sigma_cr^2*10^12); % L - Rib Spacing
    %}
end

function t = solve_skin_thickness_1(geometry, c, E, L, N, q, t_skin) % Use combined loading inequality of compressive and shear stress to solve for t % NOT TESTED

    g = geometry;
    E = E*10^9;
    
    syms t
    Rb = g.h / g.b;
    Rt = g.ts / t_skin;
    Kc = extr_zstr_k(g.h/g.b, g.ts/t_skin, g.d/g.h);
    Ks = extr_kb( L/c, 1);
    eqn = 10^-6 * t^6 * (-0.99*Ks^2*E^2)/(g.b^4) + q^2*10^6 + 10^-3*(N*Ks*E*t^3)/((1+Rt*Rb)*Kc*g.b^2) == 0;
    % eqn = 10^6*t^3 * (-0.99*Ks^2*E^2)/(g.b^4) + q^2*10^6*t + 10^3*(N*Ks*E)/((1+Rt*Rb)*Kc*g.b^2) == 0;
    S = solve(eqn, t, 'Real', true); % Solve for t
    t = eval(min(S(S>0))); % Find smallest real positive solution for t
end

function t = solve_skin_thickness_2(N, E, geometry) % based on initial buckling criteria % NOT TESTED
    
    Kc = find_Kc(Rb, Rt);
    t = ((b^2*N)/(Kc*E))^(1/3);
end

function weight = stringer_panel_weight(n, b, t, As) % NOT TESTED
    weight = b*t*(n+1) + n*As;
end

%% L varying along span
function [L_distribution, Tr_distribution, t_distribution] = varying_L_solution(geometry, E) % NOT TESTED
    g = geometry;
    F = geometry.F;
    
    wing_span = extract_dimension(0, 'wing_span');
    
    %{
    Kc = extr_zstr_k(g.h/g.b, g.ts/g.t, g.d/g.h);
    L = (3.62/Kc)^2*F^2*g.t^2*E/(10^-3*N); % Solution to Farrar Equation
    M = extract_force(L, 'BM');
    Tr = rib_thickness(geometry, c, E, bh, L, M);
    %}
    
    Tr_array = [];
    t_array = [g.t];
    L_distribution = [0];
    
    bh = extract_dimension(0, 'bh');
    c = extract_dimension(0, 'c');
    M = extract_force(0, 'BM');
    N = M/(c*bh);
    while 1 % iterate through each section in wing and assign rib spacing, rib thickness and skin thickness
        Kc = extr_zstr_k(g.h/g.b, g.ts/t_array(end), g.d/g.h);
        L = (3.62/Kc)^2*F^2*g.t^2*E/(10^-3*N); % Solution to Farrar Equation
        if L_distribution(end) + L > wing_span % check if section beyond wing_span
            break
        else
            span_location = L_distribution(end) + L;
            q = skin_shear_flow(span_location);
            M = extract_force(span_location, 'BM');
            c = extract_dimension(span_location, 'c');
            bh = extract_dimension(span_location, 'bh');
            N = M/(c*bh); % compressive Load at section*L
            
            L_distribution(end + 1) = span_location;
            Tr_array(end + 1) = rib_thickness(geometry, c, E, bh, L, M);
            t_array(end + 1) = solve_skin_thickness_1(geometry, c, E, L, N, q, t_array(end));
        end
    end
    %{
    % One last cycle for the last section past last rib
    q = skin_shear_flow(L_distribution(section));
    M = extract_force(L_distribution(section), 'BM');
    c = extract_dimension(L_distribution(section), 'c');
    bh = extract_dimension(L_distribution(section), 'bh');
    N = M/(c*bh); % compressive Load at section*L
     
    t_array(section) = solve_skin_thickness_1(N, ts/t, b/h, b, E, q);
    %}
    Tr_distribution = Tr_array;
    t_distribution = t_array;
end


%% L constant along span
function [L_distribution, Tr_distribution, t_distribution] = constant_L_solution(geometry, E) % NOT TESTED
% Calculate rib thickness based on constant L
    g = geometry;
    F = geometry.F;
    
    %load('station.mat');
    wing_span = extract_dimension(0, 'wing_span');
    
    %Tr = 5; % Rib thickness for first rib (mm) initial assumption
    M = extract_force(0, 'BM');
    bh = extract_dimension(0, 'bh');
    c = extract_dimension(0, 'c');
    N = M/(c*bh);
    
    D = bh; % rib height same as box height
    %{
    for i = 1:10 % do several iterations untill D and L converge
        Kc = extr_zstr_k(g.h/g.b, g.ts/g.t, g.d/g.h);
        sigma_cr = Kc*E*10^9*(g.t/g.b)^2; % units in N/m^2
        % sigma_cr = (Kc/3.62)* (N/g.t)/1000;
        %L = F^2*N*E*10^9/sigma_cr^2;
        L = (3.62/Kc)^2*F^2*g.t^2*70*10^9/(10^6*N); % Solution to Farrar Equation
        %L = ((4*F^2*D^2*Tr*E)/N)^(1/3)*100; % Solution for intercept of rib and stringer weight

        M = extract_force(L, 'BM');
        Tr = rib_thickness(geometry, c, E, bh, L, M);
    end
    %}
    Kc = extr_zstr_k(g.h/g.b, g.ts/g.t, g.d/g.h);
    L = (3.62/Kc)^2*F^2*g.t^2*E/(10^-3*N); % Solution to Farrar Equation
    M = extract_force(L, 'BM');
    Tr = rib_thickness(geometry, c, E, bh, L, M);
    
    Tr_array = zeros(1, floor(wing_span/L));
    t_array = zeros(1, floor(wing_span/L));
    L_distribution = (1:floor(wing_span/L))*L;
        
    for section = 1:floor(wing_span/L) % iterate through each section in wing and assign rib thickness and skin thickness
        q = skin_shear_flow(L_distribution(section));
        M = extract_force(L_distribution(section), 'BM');
        c = extract_dimension(L_distribution(section), 'c');
        bh = extract_dimension(L_distribution(section), 'bh');
        N = M/(c*bh); % compressive Load at section*L
        
        t_array(section) = solve_skin_thickness_1(geometry, c, E, L, N, q);
        D_array(section) = rib_thickness(geometry, c, E, bh, L, M);
    end
    %{
    % One last cycle for the last section past last rib
    q = skin_shear_flow(L_distribution(section));
    M = extract_force(L_distribution(section), 'BM');
    c = extract_dimension(L_distribution(section), 'c');
    bh = extract_dimension(L_distribution(section), 'bh');
    N = M/(c*bh); % compressive Load at section*L
     
    t_array(section) = solve_skin_thickness_1(N, ts/t, b/h, b, E, q);
    %}
    Tr_distribution = Tr_array;
    t_distribution = t_array;
end

function q = skin_shear_flow(L_span) % returns q (N/mm)
% shear flow is purely from torque on the top and bottom skin
    T = extract_force(L_span,'T');
	c = extract_dimension(L_span, 'c');
    b2 = extract_dimension(L_span, 'bh');
    q = T/(2*b2*c*1000);  
end

function I = second_moment_area(c, te, hc) % c:(m),te(mm),hc(mm),I(m^4)
    % hc = beam depth between panels = box height
    I = c*(te/1000)^3/12 + c*(te/1000)*(hc/2)^2;
end

function Tr = rib_thickness(geometry, c, E, hc, L, M) % Returns Tr in (mm)
% Calculates the rib thickness based on parameters
% COULD TAKE INTO ACOUNT YIELDING CRITERIA IN FUTURE
    % E in terms of GPa
    % hc is height of the wing box (hb)
    g = geometry;
    s = L;
    E = E*10^9;

    te = (g.t + g.As/g.b)/1000;
    I = second_moment_area(c, te*1000, hc);  % Second moment of area of the skin
    F_crush = (M^2*s*hc*te*c)/(2*E*I^2);
    Tr = ((F_crush*hc^2)/(3.62*c*E))^(1/3)*1000;
end