function print_factored_struct(Q_factored)
% PRINT_FACTORED_STRUCT Pretty-prints a Q_factored struct array.
%   Each entry is displayed as (1/scalar_str) * [int_vec]
%
% Usage:
%   print_factored_struct(orthonormal_set_GSP)

    fprintf('\n  { ');
    n = length(Q_factored);
    for k = 1:n
        int_vec    = Q_factored(k).int_vec;
        scalar_str = Q_factored(k).scalar_str;

        fprintf('(1/%s) * [', scalar_str);
        for r = 1:length(int_vec)
            if r == 1
                fprintf('%g', int_vec(r));
            else
                fprintf(', %g', int_vec(r));
            end
        end
        fprintf(']');

        if k < n
            fprintf(',\n    ');
        end
    end
    fprintf(' }\n');
end
