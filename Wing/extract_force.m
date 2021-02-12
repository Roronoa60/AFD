% This function returns either the shear force (SF), bending moment (BM) or torque (T)  at a specified span

function force = extract_force(span_location, force_type)
    load('station.mat');
    % load('force_distributions.mat');
    span_end = SpanMesh(end);
    
    if 0 < span_location <= span_end
        span_index = find((span_location-0.005 < SpanMesh) & (SpanMesh < span_location+0.005));
        if force_type == 'SF'
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