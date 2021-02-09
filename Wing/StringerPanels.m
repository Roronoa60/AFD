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
function  [n, b, t2, ts, h, d] = root_geometry(M, c, bh, bh_min, E, q) %(n, t2, ts, h, M, c, b2, E)
    % M - root bending moment (Nm)
    % c - root box width (m)
    % bh - root box height (m)
    % bh_min - min box height
    % E - Young Modulus of panel and stringer material (Pa)
    % q - Shear Flow

    N = M/(c*bh); % N - Compressive Load of unit length (N/m)
    
    weights = {'n', 'b', 't2', 't2', 'ts', 'h', 'd', 'F', 'weight'};
    min_weight = {'n', 'b', 't2', 't2', 'ts', 'h', 'd', 'F', 'weight'};
    
    % n - No. of stringers
    % Rt = ts/t = stringer web thickness/skin thickness
    % Rb = d/b = stringer height/stringer pitch  !!!!! d string height/flange width
    % d_h = d/h
    for n = 1:16
        for Rt = [0.5,0.6,0.7,0.9,1,1.25,1.5,2]
            for d_h = [0.3,0.4,0.5]
                for h_b = 1 %!!!!
                    % Check h < bh_min/2
                    As = (h + 2*d)*ts; %As stringer area
                    % F = FARRAR(Rt, As/(b*t)) IS t=ts or t!!!!!!
                    % t = solve_skin_thickness(N, d_h, Rt, Rb, h_b, E)
                    % update min_weight
                    % weights(end+1,:) = {'n', 'b', 't2', 't2', 'ts', 'h', 'd', 'F', 'weight'}
        
                end
            end
        end
    end
                    
        



    
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
    
end

%% 3.2.6
