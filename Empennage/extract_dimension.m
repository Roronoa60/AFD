% function to extract wing box dimensions at a given span location

function dimension = extract_dimension(span_location, dimension_type, VT_or_HT)
    switch VT_or_HT
        case 'VT'
            load('station_Vtail.mat');
        case 'HT'
            load('station_htp.mat');
    end
    span_end = SpanMesh(end);
    if 0 < span_location <= span_end
        span_index = min(find((span_location-0.005 < SpanMesh) & (SpanMesh <= span_location+0.0051)));
        switch dimension_type
            case 'c'
                dimension = c(span_index);
            case 'bh'
                dimension = b2(span_index);
            case 'chord'
                dimension = Chord_len(span_index);
            case 'wing_span'
                dimension = span_end;
            otherwise
                error('wing box dimension not recognized');
        end
    else
        error('span_location needs to be within span limit');
    end  
end