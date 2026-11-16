# File formats: ALCDEF light-curve blocks, CdR (CourbRot) files; MPC band corrections; light-time correction.

const AU_LIGHT_DAYS = 499.004783836 / 86400     # light time for 1 au [d]

"""
    light_time_correct(jd, Δ) -> jd − Δ τ_au

Times corrected for the light travel time from the object (geocentric distance `Δ` [au]) to the observer.
"""
light_time_correct(jd, Δ) = jd .- Δ .* AU_LIGHT_DAYS

"""
    MPC_BAND_TO_V

Approximate corrections to the V band used by the Minor Planet Center for reported magnitudes, V = m + Δ(band)
(average asteroid colours; e.g. "R" → +0.4, "G" (Gaia) → +0.28, "w" (Pan-STARRS) → −0.13). A blank band is "".
For objects with known colours use a proper transformation (e.g. [`ps1_to_V`](@ref)) or fit band offsets
([`fit_phase_curve`](@ref)).
"""
const MPC_BAND_TO_V = Dict{String,Float64}(
    "" => -0.8, "U" => -1.3, "B" => -0.8, "g" => -0.35, "V" => 0.0, "r" => 0.14, "R" => 0.4, "C" => 0.4, "W" => 0.4,
    "i" => 0.32, "z" => 0.26, "I" => 0.8, "J" => 1.2, "w" => -0.13, "y" => 0.32, "L" => 0.2, "H" => 1.4, "K" => 1.7,
    "Y" => 0.7, "G" => 0.28, "v" => 0.0, "c" => -0.05, "o" => 0.33, "u" => 2.5)

"""
    mpc_to_V(mag, band) -> V

V magnitude from an MPC-reported magnitude with the band corrections [`MPC_BAND_TO_V`](@ref). Pan-STARRS bands
reported as "Pg", "Pr", "Pw", … are mapped to "g", "r", "w", ….
"""
function mpc_to_V(mag::Real, band::AbstractString)
    b = length(band) == 2 && band[1] == 'P' ? string(band[2]) : String(band)
    haskey(MPC_BAND_TO_V, b) || throw(ArgumentError("no MPC V-band correction for band \"$band\""))
    return mag + MPC_BAND_TO_V[b]
end

# ------------------------------------------------------------------------------------------------ ALCDEF

"""
    ALCDEFBlock

One ALCDEF light-curve block (ALCDEF standard v2.3, Warner et al.): `meta`, the keyword/value pairs of the metadata
section in file order (keywords may repeat, e.g. `EQUINOX`), and the data columns `jd` (UT Julian date of
mid-exposure), `mag` and `err` (`NaN` if not given). `block["KEYWORD"]` returns the first value of a keyword or
`nothing`; `block["KEYWORD", default]` with a default.
"""
struct ALCDEFBlock
    meta::Vector{Pair{String,String}}
    jd::Vector{Float64}
    mag::Vector{Float64}
    err::Vector{Float64}
end
function Base.getindex(b::ALCDEFBlock, key::AbstractString)
    i = findfirst(p -> p.first == uppercase(key), b.meta)
    return i === nothing ? nothing : b.meta[i].second
end
Base.getindex(b::ALCDEFBlock, key::AbstractString, default) = (v = b[key]; v === nothing || isempty(v) ? default : v)
Base.length(b::ALCDEFBlock) = length(b.jd)
function Base.show(io::IO, b::ALCDEFBlock)
    print(io, "ALCDEFBlock(", b["OBJECTNUMBER", "?"], " ", b["OBJECTNAME", ""], ", ", b["SESSIONDATE", "?"], " ",
          b["SESSIONTIME", ""], ", ", length(b), " points, ", b["FILTER", "?"], "/", b["MAGBAND", "?"], ", ",
          b["CONTACTNAME", ""], ")")
end

const _ALCDEF_DELIM = Dict("PIPE" => '|', "SEMIC" => ';', "SEMICOLON" => ';', "TAB" => '\t', "COMMA" => ',', "SPACE" => ' ')

"""
    read_alcdef(file_or_io) -> Vector{ALCDEFBlock}

Read all light-curve blocks of an ALCDEF file (STARTMETADATA … ENDMETADATA, DATA=JD|MAG[|MAGERR[|AIRMASS]] …
ENDDATA). Blocks without `ENDDATA` (as in the ALCDEF bulk archive) end at the next `STARTMETADATA` or at the end of
the file. The data delimiter is taken from the `DELIMITER` keyword (PIPE, SEMIC, TAB; COMMA and SPACE of the simple
ALCDEF are accepted); CRLF line ends and blank lines are tolerated.
"""
read_alcdef(file::AbstractString) = open(read_alcdef, file)
function read_alcdef(io::IO)
    blocks = ALCDEFBlock[]
    meta = Pair{String,String}[]; jd = Float64[]; mag = Float64[]; err = Float64[]
    state = :outside
    delim = '|'
    for raw in eachline(io)
        line = strip(raw)
        isempty(line) && continue
        if line == "STARTMETADATA"
            state == :data && push!(blocks, ALCDEFBlock(meta, jd, mag, err))     # block without ENDDATA
            meta = Pair{String,String}[]; jd = Float64[]; mag = Float64[]; err = Float64[]
            state = :meta; delim = '|'
        elseif line == "ENDMETADATA"
            state = :data
        elseif line == "ENDDATA"
            state == :data && push!(blocks, ALCDEFBlock(meta, jd, mag, err))
            state = :outside
        elseif state == :meta
            i = findfirst('=', line)
            i === nothing && continue
            k = uppercase(strip(line[1:prevind(line, i)])); v = strip(line[nextind(line, i):end])
            push!(meta, k => String(v))
            k == "DELIMITER" && (delim = get(_ALCDEF_DELIM, uppercase(v), '|'))
        elseif state == :data && startswith(line, "DATA=")
            f = split(line[6:end], delim; keepempty=true)
            length(f) >= 2 || continue
            push!(jd, parse(Float64, strip(f[1]))); push!(mag, parse(Float64, strip(f[2])))
            push!(err, length(f) >= 3 && !isempty(strip(f[3])) ? parse(Float64, strip(f[3])) : NaN)
        end
    end
    state == :data && push!(blocks, ALCDEFBlock(meta, jd, mag, err))             # last block without ENDDATA
    return blocks
end

"""
    write_alcdef(file_or_io, blocks)

Write ALCDEF light-curve blocks (metadata in the stored order, data with the PIPE delimiter; `DELIMITER=PIPE` is
enforced). Errors that are `NaN` are omitted.
"""
write_alcdef(file::AbstractString, blocks) = open(io -> write_alcdef(io, blocks), file, "w")
function write_alcdef(io::IO, blocks::AbstractVector{ALCDEFBlock})
    for b in blocks
        println(io, "STARTMETADATA")
        hasdelim = false
        for (k, v) in b.meta
            if k == "DELIMITER"
                println(io, "DELIMITER=PIPE"); hasdelim = true
            else
                println(io, k, "=", v)
            end
        end
        hasdelim || println(io, "DELIMITER=PIPE")
        println(io, "ENDMETADATA")
        for i in eachindex(b.jd)
            if isnan(b.err[i])
                @printf(io, "DATA=%.6f|%+.3f\n", b.jd[i], b.mag[i])
            else
                @printf(io, "DATA=%.6f|%+.3f|%.3f\n", b.jd[i], b.mag[i], b.err[i])
            end
        end
        println(io, "ENDDATA")
    end
end
write_alcdef(io::IO, b::ALCDEFBlock) = write_alcdef(io, [b])

"""
    alcdef_table(blocks) -> NamedTuple

All data points of a set of ALCDEF blocks as columns `(jd, mag, err, session)`, where `session` is the block index
(use it as `groups` in [`falc`](@ref) / [`fit_period`](@ref) for differential magnitudes with one zero point per
session), and `labels`, one label per block ("SESSIONDATE CONTACTNAME").
"""
function alcdef_table(blocks::AbstractVector{ALCDEFBlock})
    jd = reduce(vcat, [b.jd for b in blocks]); mag = reduce(vcat, [b.mag for b in blocks])
    err = reduce(vcat, [b.err for b in blocks])
    session = reduce(vcat, [fill(i, length(b)) for (i, b) in enumerate(blocks)])
    labels = [string(b["SESSIONDATE", "?"], " ", b["CONTACTNAME", ""]) for b in blocks]
    return (jd=jd, mag=mag, err=err, session=session, labels=labels)
end

# ------------------------------------------------------------------------------------------------ CdR (CourbRot)

"""
    read_cdr(file) -> NamedTuple

Read a light-curve file in the CdR (CourbRot) format of R. Behrend (header lines with three-letter keywords `FMT`,
`NOM`, `MES`, `POS`, `TEL`, …, then whitespace-separated data; lines starting with `#` or `;` are comments, e.g.
rejected points). The `FMT` string gives the column meaning; `D` (date, JD) and `V` (target magnitude) are required,
`R` columns are magnitudes of comparison stars. `POS p e` gives the time stamp position (−1 start, 0 middle, +1 end of
the exposure) and the exposure time `e` [s]; returned times are mid-exposure. Returns
`(jd, mag, dmag, comps, object, observer, header::Dict, columns::Vector{Vector{String}})`, where `comps` is the matrix
of comparison-star magnitudes (one column per `R`) and `dmag` the differential magnitude V − mean(R) (equal to `mag`
without `R` columns), as in the CdR / AsPyLib light-curve tools.
"""
function read_cdr(file::AbstractString)
    header = Dict{String,String}()
    rows = Vector{Vector{String}}()
    for raw in eachline(file)
        line = strip(raw)
        (isempty(line) || line[1] in ('#', ';')) && continue
        if occursin(r"^[A-Za-z]{3}(\s|$)", line)
            header[uppercase(line[1:3])] = strip(line[4:end])
        else
            push!(rows, String.(split(line)))
        end
    end
    fmt = replace(get(header, "FMT", "DV"), " " => "")
    iD = findfirst(==('D'), fmt); iV = findfirst(==('V'), fmt)
    (iD === nothing || iV === nothing) && throw(ArgumentError("$file: FMT must contain D and V"))
    jd = [parse(Float64, r[iD]) for r in rows]; mag = [parse(Float64, r[iV]) for r in rows]
    iR = findall(==('R'), fmt)
    comps = [parse(Float64, r[j]) for r in rows, j in iR]
    dmag = isempty(iR) ? copy(mag) : mag .- vec(mean(comps; dims=2))
    pos = split(get(header, "POS", "0 0"))
    p = parse(Float64, pos[1]); e = length(pos) > 1 ? parse(Float64, pos[2]) : 0.0
    jd .-= p * e / 2 / 86400
    return (jd=jd, mag=mag, dmag=dmag, comps=comps, object=get(header, "NOM", ""), observer=get(header, "MES", ""),
            header=header, columns=rows)
end
