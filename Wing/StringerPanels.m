clear
clc

[t, sigma] = simple_panel_buckling(23286.8, 1, 0.2, 70E9);

%stringer_panel_buckling_thickness(10, 1, 1, 50, 23286.8, 1, 0.2, 70E9);


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
    weights = {'n', 'b', 't', 'ts', 'h', 'd', 'F', 'As', 'weight'};
    min_weight = {'n', 'b', 't', 'ts', 'h', 'd', 'F', 'As', 0};
    max_F = {'n', 'b', 't', 'ts', 'h', 'd', 0, 'As', 'weight'};
    
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
                        if min_weight(end) < weight_save
                            min_weight = {n, b, t, ts, h, d, F, As, weight_save};
                        end
                        % If new geometry improves F
                        if max_F(7) < F | (max_F(7) == F & max_F(end) > weight_save) % if F matches, pick lighter geometry
                            max_F = {n, b, t, ts, h, d, F, As, weight_save};
                        end
                        % weights(end+1,:) = {n, b, t, ts, h, d, F, weight_save}
                    else
                        % print out some message that h is too large?
                    end
                end
            end
        end
    end
    min_weight_geometry = min_weight;
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

function t = solve_skin_thickness_1(N, Rt, Rb, b, E, q) % Use combined loading inequality of compressive and shear stress to solve for t
    %{
    syms t

    eqn = (t^4 - 2*t)*(t-2) == 0;
    S = solve(eqn,t,'Real',true)
    S = min(S(S>0))
    %}
    syms t
    
    Ks = find_Ks(Rb, Rt);
    Kc = find_Kc(Rb, Rt);
    eqn = t^6 * (-0.99*Ks^2*E^2)/(b^4) + q^2 + (N*Ks*E*t^3)/((1+Rt*Rb)*Kc*b^2) == 0;
    
    S = solve(eqn, t, 'Real', true); % Solve for t
    t = min(S(S>0)); % Find smallest real positive solution for t
end

function t = solve_skin_thickness_2(b, N, E, Rb, Rt) % based on initial buckling criteria
    Kc = find_Kc(Rb, Rt);
    t = ((b^2*N)/(Kc*E))^(1/3);
end

function weight = stringer_panel_weight(n, b, t, As)
    weight = b*t*(n+1) + n*As;
end

%% L varying along span
function [L, D] = varying_L_solution(h, b, ts, t, E, N)
    % Obtain optimum L using Farrer's Equation
    Kc = get_Kc(h/b, ts/t);
    sigma_cr = Kc*E*(t/b)^2;
    L = F^2*N*E/sigma_cr; % From rearanging FARRAR's equation
end
% Calculate L for each station
% Calcuate rib thickness based on rib spacing
% Decide whether calculated thickness B is manufacturable
% Caculate weight of entire configuration

%% L constant along span
function [L, D_distribution] = constant_L_solution(F, E, bh, ts, t, As, b, h, c wing_span)
    % Calculate rib thickness based on constant L
    D = 5; % (mm) initial assumption
    M = bending_moment_at('0_span');
    N = M/(c*bh);
    
    for i = 1:10 % do several iterations untill D and L converge
        L = ((4*F^2*D^2*Tr*E)/N)^(1/3); % Solution for intercept of rib and stringer weight

        M = bending_moment_at_section(L); % !!!!!!!!!!!!!
        D = rib_thickness(t, As, b, c, te, hb, M, L, hc, E);
    end
    D_array = zeros(1, round(wing_span/L));
    t_array = zeros(1, round(wing_span/L)+1);
    D_array(1) = D;
    t_array(1) = t;
    section = 2;
    while section*L < wing_span % iterate through each section in wing and assign rib thickness and skin thickness
        q = shear_flow_at(section*L);
        M = compressive_load_at(section*L);
        c = c_at(section*L);
        bh = bh_at(section*L);
        N = M/(c*bh); % compressive Load at section*L
        
        t_array(section) = solve_skin_thickness_1(N, ts/t, b/h, b, E, q);
        D_array(section) = rib_thickness(t, As, b, c, hb, M, L, hc, E);
        section = section + 1;
    end
    q = shear_flow_at(section*L);
    M = compressive_load_at(section*L);
    c = c_at(section*L);
    bh = bh_at(section*L);
    N = M/(c*bh); % compressive Load at section*L
        
    t_array(section) = solve_skin_thickness_1(N, ts/t, b/h, b, E, q);
        
end


function I = second_moment_area(c, te, hc)
    % hc = beam depth between panels = box height
    I = c*(te/1000)^3/12 + c*(te/1000)*(hc/2)^2;
end

function D = rib_thickness(t, As, b, c, hb, M, s, hc, E)
        te = t + As/b;
        I = second_moment_area(c, te, hb);  % Second moment of area of the skin
        F_crush = (M^2*s*hc*te*c)/(2*E*I^2);
        D = ((F_crush*hc^2)/(3.62*c*E))^(1/3);
end