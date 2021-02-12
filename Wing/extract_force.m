% This function returns either the shear force (SF), bending moment (BM) or torque (T)  at a specified span

function force = extract_force(span_location, force_type)
    load('station.mat');
    load('force_distributions.mat');
    span_end = station.SpanMesh(end);
    
    if 0 < span_location <= span_end
        span_index = find(station.SpanMesh == round(span_location,2));
        if force_type == 'SF'
            force = force.SF(span_index);
        elseif force_type == 'BM'
            force = force.BM(span_index);
        elseif force_type == 'T'
            force = force.T(span_index);
        else
            error('force_type not recognized')
        end
    else
        error('span_location needs to be within span limit');
    end       
    
end