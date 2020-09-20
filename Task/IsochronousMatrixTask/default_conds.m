function [selections] = default_conds(default_selection)

switch default_selection
    case '2.4Hz-isochronous-matrix'
        iso_mat_24
    otherwise
        select_conditions_GUI
end


end