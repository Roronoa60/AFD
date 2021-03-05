% This function returns either the shear force (SF), bending moment (BM) or torque (T)  at a specified span

function force = extract_force(span_location, force_type, loading_direction)
    switch loading_direction
        case 'top'
            load('station_ult.mat');
        case 'bottom'
            load('station_negative');
    end
    
    % load('force_distributions.mat');
    span_end = SpanMesh(end);
    
    if 0 < span_location <= span_end
        span_index = find((span_location-0.005 < SpanMesh) & (SpanMesh < span_location+0.005));
        if force_type == 'SF'
            if loading_direction == 'top'
                force_ult = SF_tot(span_index);
                load('station_land.mat')
                force_land = SF_tot(span_index);
                if force_ult > force_land
                    force = force_ult;
                else
                    force = force_land;
                end
            end
            force = SF_tot(span_index);
        elseif force_type == 'BM'
            force = BM_tot(span_index);
        elseif force_type == 'T'
            force = T(span_index);
        else
            error('force_type not recognized')
        end
    else
        error('span_location needs to be within span limit');
    end       
    
end