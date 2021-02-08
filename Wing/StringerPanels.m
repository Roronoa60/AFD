clear
clc

simple_panel_buckling(23286.8, 1, 0.2, 70E9);

stringer_panel_buckling_thickness(10, 1, 1, 50, 23286.8, 1, 0.2, 70E9)

function  [t2,sigma] = simple_panel_buckling(M, c, b2, E)
    % M - Bending Moment at section (Nm)
    % c - width (of box or panel section) (m)
    % b2 - Web height (m) 
    % E - Young Modulus of panel (Pa)
    
    N = M/(c*b2); % N - Compressive Load of unit length (N/m)
    t2 = (N/(3.62*E)*c^2)^(1/3)*1000; % t2 - Skin Thickness (mm)
    sigma = N/t2/1000; % sigma - Buckling stress (N/mm^2)
end

function  output = stringer_panel_buckling_thickness(n, t2, ts, h, M, c, b2, E)
    % n - No. of stringers
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
    sigma_ratio = 1.38 % sigma_ratio = sigma_cr/sigma_0 = Catchpole_fun(ts/t2,As/(b*t2))   READ FROM CATCHPOLE DIAGRAM
     
    % FARRAR efficiency factor
    F = 0.82; % F = FARRAR_Factor(ts/t2, As/(b*t2))
    % sigma_f = F/sqrt(L/N/E)/10^6 % sigma_f - Mean stress by skin and stringer
    sigma_cr = sigma_0*sigma_ratio % sigma_cr - critical buckling stress
    
    % Equate critical initial buckling stress with critical flexural
    % buckling and rearange for L:
    L = F^2*N*E/(sigma_cr^2*10^12) % L - Rib Spacing
    
end
