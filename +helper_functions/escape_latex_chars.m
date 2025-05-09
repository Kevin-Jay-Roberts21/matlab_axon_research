% Helper function to catch unexpected characters in data names
% Kevin Roberts
% May 2025

function escaped = escape_latex_chars(str)
    special_chars = {'\', '_', '^', '{', '}', '$', '%', '#', '&', '~'};
    replacements = {'\\textbackslash{}', '\_', '\^{}', '\{', '\}', '\$', '\%', '\#', '\&', '\\textasciitilde{}'};

    escaped = str;
    for k = 1:length(special_chars)
        escaped = strrep(escaped, special_chars{k}, replacements{k});
    end
end