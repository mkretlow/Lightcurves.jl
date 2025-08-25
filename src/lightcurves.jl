"""
## Purpose
    Functions for asteroid light curves here.

## Functions
    ellipsoid_full_amp : 
"""

export ellipsoid_full_amp


function ellipsoid_full_amp(a,b,c,aa)

    a = (1 - cos(aa)^2) * (a^2/b^2 - 1)
    b = 1 + cos(aa)^2 * (a^2/b^2 - 1)
    
    D = a / b
    
    return 1.25 * log10(1 + D)
end


# Fit a given period estimate/guess
function fit_period(xdata, ydata, weights, nb_harmon)

    nb_harmon = 4 #6

    an = ones(nb_harmon)
    bn = ones(nb_harmon)

    p = vcat(an,bn)

    function model(x, p, w, nb_harmon)
        result = p[1]

        for i=1:nb_harmon
            result = result + p[i]*cos(i*w*x) + p[i+nb_harmon]*sin(i*w*x)
        end

        return result
    end

    p = vcat(p,w)

    fit = curve_fit(model, xdata, ydata, p)

end

