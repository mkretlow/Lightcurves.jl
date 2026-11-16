# FITS image reading for Lightcurves.readfits (loaded with `using FITSIO`).
module LightcurvesFITSIOExt

using Lightcurves
using FITSIO

function _readfits(file::AbstractString, hdu)
    FITS(file) do io
        k = hdu
        if k === nothing
            k = findfirst(h -> h isa ImageHDU && length(size(h)) == 2, collect(io))
            k === nothing && error("$file: no HDU with a 2-D image")
        end
        h = io[k]
        img = read(h)
        ndims(img) == 2 || error("$file, HDU $k: image has $(ndims(img)) dimensions, expected 2")
        return Float64.(img), read_header(h)
    end
end

__init__() = (Lightcurves._FITS_READER[] = _readfits)

end
