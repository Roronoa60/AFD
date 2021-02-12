% function to extract wing box dimensions at a given span location

function dimension = extract_dimension(span_location, dimension_type)
    load('station.mat');
    span_end = station.SpanMesh(end);
    
    if 0 < span_location <= span_end
        span_index = find(station.SpanMesh == round(span_location,2));
        switch dimension_type
            case 'c'
                dimension = station.Chord(span_index);
            case 'bh'
                dimension = station.bh(span_index);
            otherwise
                error('wing box dimension not recognized');
        end
    else
        error('span_location needs to be within span limit');
    end  
end