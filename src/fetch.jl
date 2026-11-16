# Data access: MPC observations (photometry), JPL Horizons geometry, JPL SBDB, ALCDEF light curves, LCDB.
# All requests send only the object designation; no user information. Large downloads are avoided (single members of
# the LCDB zip are read with HTTP range requests).

const USER_AGENT = "Lightcurves.jl (Julia)"

"""
    cachedir() -> String

Directory for downloaded files: `ENV["LIGHTCURVES_CACHE"]` or `~/.cache/Lightcurves.jl`.
"""
cachedir() = get(ENV, "LIGHTCURVES_CACHE", joinpath(homedir(), ".cache", "Lightcurves.jl"))

"""
    datadir() -> String

Directory of the local databases (ALCDEF archive, LCDB release and the files extracted from it):
`ENV["LIGHTCURVES_DATA"]` or the `data/` folder of the package. Other downloads are cached in [`cachedir`](@ref).
"""
datadir() = get(ENV, "LIGHTCURVES_DATA", joinpath(dirname(@__DIR__), "data"))

function _http(url; method="GET", body=nothing, headers=Pair{String,String}[], timeout=120)
    out = IOBuffer()
    hdr = vcat(["User-Agent" => USER_AGENT], headers)
    resp = Downloads.request(url; method=method, output=out, headers=hdr, timeout=timeout,
                             input=body === nothing ? nothing : IOBuffer(body), throw=false)
    resp isa Downloads.Response || throw(ErrorException("request to $url failed: $resp"))
    200 <= resp.status < 300 || resp.status == 206 ||
        throw(ErrorException("HTTP $(resp.status) from $url"))
    return take!(out), resp
end
_get(url; kw...) = String(first(_http(url; kw...)))

_urlencode(s) = join(c in ('-', '.', '_', '~') || isletter(c) && isascii(c) || isdigit(c) ? string(c) :
                     c == ' ' ? "+" : join("%" * uppercase(string(b; base=16, pad=2)) for b in codeunits(string(c))) for c in s)
_form(pairs) = join([_urlencode(k) * "=" * _urlencode(v) for (k, v) in pairs], "&")
_post(url, pairs) = String(first(_http(url; method="POST", body=_form(pairs),
                                       headers=["Content-Type" => "application/x-www-form-urlencoded"])))

# ------------------------------------------------------------------------------------------------ MPC

"""
    fetch_mpc_observations(desig) -> Vector{Dict}

All observations of an object from the MPC observations API (`get-obs`, ADES format), e.g. `"612533"` or `"2002 XV93"`.
Each element is a Dict of ADES fields (`obstime`, `stn`, `mag`, `band`, `rmsmag`, `photcat`, `ra`, `dec`, …).
"""
function fetch_mpc_observations(desig::AbstractString)
    body = JSON.json(Dict("desigs" => [desig], "output_format" => ["ADES_DF"]))
    txt = String(first(_http("https://data.minorplanetcenter.net/api/get-obs"; method="GET", body=body,
                             headers=["Content-Type" => "application/json"])))
    d = JSON.parse(txt)
    isempty(d) && return Dict{String,Any}[]
    obs = d[1]["ADES_DF"]
    return obs === nothing ? Dict{String,Any}[] : [Dict{String,Any}(String(k) => v for (k, v) in o) for o in obs]
end

"""
    fetch_mpc_photometry(desig) -> NamedTuple

The observations with a magnitude: `(jd, mag, band, rmsmag, stn, photcat)` (`rmsmag` NaN if not reported; JD UTC).
"""
function fetch_mpc_photometry(desig::AbstractString)
    obs = filter(o -> get(o, "mag", nothing) !== nothing && !isempty(string(o["mag"])), fetch_mpc_observations(desig))
    val(o, k) = (v = get(o, k, nothing); v === nothing ? "" : string(v))
    jd = [_iso_to_jd(val(o, "obstime")) for o in obs]
    return (jd=jd, mag=[parse(Float64, val(o, "mag")) for o in obs], band=[val(o, "band") for o in obs],
            rmsmag=[isempty(val(o, "rmsmag")) ? NaN : parse(Float64, val(o, "rmsmag")) for o in obs],
            stn=[val(o, "stn") for o in obs], photcat=[val(o, "photcat") for o in obs])
end

function _iso_to_jd(s::AbstractString)
    s = replace(s, "Z" => "")
    m = match(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2}(?:\.\d+)?)", s)
    m === nothing && throw(ArgumentError("cannot parse time $s"))
    y, mo, d, h, mi = parse.(Int, m.captures[1:5]); sec = parse(Float64, m.captures[6])
    return Dates.datetime2julian(DateTime(y, mo, d, h, mi)) + sec / 86400
end

# ------------------------------------------------------------------------------------------------ JPL

"""
    fetch_horizons_geometry(desig; start, stop, step="1d", center="500") -> NamedTuple

Observing geometry from JPL Horizons (observer table; `center="500"` geocentric, or an MPC station code as
`"code@399"`): `(jd, r, delta, alpha, vmag)` with heliocentric and observer distance [au], solar phase angle [deg] and
the predicted magnitude. `start`, `stop` as "YYYY-MM-DD"; `step` e.g. "1d", "6h". The object is looked up by number or
designation (`desig;` small-body search).
"""
function fetch_horizons_geometry(desig::AbstractString; start::AbstractString, stop::AbstractString, step::AbstractString="1d",
                                 center::AbstractString="500")
    q = ["format" => "text", "COMMAND" => "'$(desig);'", "OBJ_DATA" => "'NO'", "MAKE_EPHEM" => "'YES'",
         "EPHEM_TYPE" => "'OBSERVER'", "CENTER" => "'$center'", "START_TIME" => "'$start'", "STOP_TIME" => "'$stop'",
         "STEP_SIZE" => "'$step'", "QUANTITIES" => "'9,19,20,43'", "CSV_FORMAT" => "'YES'"]
    txt = _get("https://ssd.jpl.nasa.gov/api/horizons.api?" * _form(q))
    return parse_horizons_geometry(txt)
end

"""
    parse_horizons_geometry(text) -> NamedTuple

Parse the CSV observer table of a Horizons request with QUANTITIES 9, 19, 20, 43 (see
[`fetch_horizons_geometry`](@ref)).
"""
function parse_horizons_geometry(txt::AbstractString)
    L = split(txt, '\n')
    i0 = findfirst(l -> strip(l) == "\$\$SOE", L); i1 = findfirst(l -> strip(l) == "\$\$EOE", L)
    (i0 === nothing || i1 === nothing) && throw(ErrorException("Horizons: no ephemeris in the response:\n" * first(txt, 600)))
    mon = Dict(m => i for (i, m) in enumerate(("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))
    jd = Float64[]; r = Float64[]; Δ = Float64[]; α = Float64[]; V = Float64[]
    for l in L[i0+1:i1-1]
        p = strip.(split(l, ','))
        m = match(r"(\d{4})-(\w{3})-(\d{2}) (\d{2}):(\d{2})(?::(\d{2}(?:\.\d+)?))?", p[1])
        m === nothing && continue
        sec = m.captures[6] === nothing ? 0.0 : parse(Float64, m.captures[6])
        push!(jd, Dates.datetime2julian(DateTime(parse(Int, m.captures[1]), mon[m.captures[2]], parse(Int, m.captures[3]),
                                                 parse(Int, m.captures[4]), parse(Int, m.captures[5]))) + sec / 86400)
        num(x) = (y = tryparse(Float64, x); y === nothing ? NaN : y)
        push!(V, num(p[4])); push!(r, num(p[6])); push!(Δ, num(p[8])); push!(α, num(p[10]))
    end
    return (jd=jd, r=r, delta=Δ, alpha=α, vmag=V)
end

"""
    observing_geometry(desig, jd; step="1d", center="500") -> NamedTuple

Heliocentric distance, observer distance and phase angle at the Julian dates `jd`, interpolated linearly in a
Horizons table covering them: `(r, delta, alpha)`. For distant objects a daily step is plenty.
"""
function observing_geometry(desig::AbstractString, jd::AbstractVector{<:Real}; step="1d", center="500")
    d(j) = Dates.format(Dates.julian2datetime(j), "yyyy-mm-dd")
    g = fetch_horizons_geometry(desig; start=d(minimum(jd) - 2), stop=d(maximum(jd) + 2), step=step, center=center)
    return (r=_interp(g.jd, g.r, jd), delta=_interp(g.jd, g.delta, jd), alpha=_interp(g.jd, g.alpha, jd))
end
function _interp(x, y, xi)
    map(xi) do v
        k = searchsortedlast(x, v)
        (k < 1 || k >= length(x)) && throw(ArgumentError("time outside the table"))
        f = (v - x[k]) / (x[k+1] - x[k]); y[k] + f * (y[k+1] - y[k])
    end
end

"""
    fetch_sbdb(desig) -> Dict

The JPL Small-Body Database record (`sbdb.api` with physical parameters): returns a Dict with `"fullname"`, `"orbit"`
(elements) and `"phys"`, a Dict of physical parameters name => (value, sigma, reference) (e.g. "H", "diameter",
"albedo", "rot_per").
"""
function fetch_sbdb(desig::AbstractString)
    d = JSON.parse(_get("https://ssd-api.jpl.nasa.gov/sbdb.api?" * _form(["sstr" => desig, "phys-par" => "1"])))
    haskey(d, "object") || throw(ErrorException("SBDB: object $desig not found"))
    phys = Dict{String,Any}()
    for p in get(d, "phys_par", [])
        phys[String(p["name"])] = (value=p["value"], sigma=get(p, "sigma", nothing), ref=get(p, "ref", nothing))
    end
    return Dict{String,Any}("fullname" => d["object"]["fullname"], "orbit" => get(d, "orbit", nothing), "phys" => phys)
end

# ------------------------------------------------------------------------------------------------ ALCDEF

const ALCDEF_URL = "https://www.alcdef.org/php/alcdef_GenerateALCDEFPage.php"

_html_text(s) = strip(replace(replace(replace(s, r"<[^>]+>" => ""), "&nbsp;" => " "), "&amp;" => "&"))

function _alcdef_objectpage(object::AbstractString)
    p1 = _post(ALCDEF_URL, ["findValue" => object, "WholeWord" => "TRUE", "findObjects" => "Submit"])
    cands = [m.captures[1] for m in eachmatch(r"name='selectedObject' value=\"([^\"]*)\"", p1)]
    isempty(cands) && throw(ErrorException("ALCDEF: no object found for \"$object\""))
    sel = something(findfirst(c -> lowercase(c) == lowercase(object), cands), 1)
    return _post(ALCDEF_URL, ["selectedObject" => cands[sel], "doSearch2" => "Display"]), cands
end

"""
    fetch_alcdef(object; ids=nothing) -> Vector{ALCDEFBlock}

Download the ALCDEF light curves of an object (number or name, e.g. `"1486"`) from alcdef.org, through the same
search–select–download steps as the web interface (no bulk archive). `ids` restricts the download to some light-curve
block IDs (see [`alcdef_sessions`](@ref)). The files contain the observers' contact information: do not redistribute
them unchanged.
"""
function fetch_alcdef(object::AbstractString; ids=nothing)
    page, _ = _alcdef_objectpage(object)
    allids = [m.captures[1] for m in eachmatch(r"name=lightcurves\[\] value='(\d+)'", page)]
    isempty(allids) && return ALCDEFBlock[]
    use = ids === nothing ? allids : string.(ids)
    num = match(r"name='SearchNumber' value='([^']*)'", page); nam = match(r"name='SearchName' value=\"([^\"]*)\"", page)
    pairs = Pair{String,String}["SearchNumber" => (num === nothing ? "" : num.captures[1]),
                                "SearchName" => (nam === nothing ? "" : nam.captures[1]), "saveLC" => "Download"]
    append!(pairs, ["lightcurves[]" => i for i in use])
    txt = _post(ALCDEF_URL, pairs)
    return read_alcdef(IOBuffer(txt))
end

"""
    alcdef_sessions(object) -> Vector{NamedTuple}

The list of ALCDEF light-curve blocks of an object without downloading the data: `(id, date, filter, band, phase,
pabl, pabb, contact)`.
"""
function alcdef_sessions(object::AbstractString)
    page, _ = _alcdef_objectpage(object)
    ids = [m.captures[1] for m in eachmatch(r"name=lightcurves\[\] value='(\d+)'", page)]
    rows = NamedTuple[]
    for id in ids
        i = findfirst("value='$id'", page)
        i === nothing && continue
        seg = page[last(i):min(end, last(i) + 3000)]
        cells = [_html_text(m.captures[1]) for m in eachmatch(r"<td[^>]*>(.*?)</td>"s, seg)]
        c = filter(!isempty, cells)
        length(c) >= 8 || continue
        push!(rows, (id=parse(Int, id), date=c[1], filter=c[2], band=c[3], phase=tryparse(Float64, c[4]),
                     pabl=tryparse(Float64, c[5]), pabb=tryparse(Float64, c[6]), contact=c[7]))
    end
    return rows
end

# ------------------------------------------------------------------------------------------------ LCDB

"""
    lcdb_entry(object) -> Dict{String,String}

The LCDB summary entry of one object (Warner et al., Asteroid Lightcurve Database), as shown on its ALCDEF page:
"Number", "Name", "Desig", "Class", "Diam", "H", "Albedo", "Period", "AMin", "AMax", "U" (quality code), … . No bulk
download. Empty if the object is not in the LCDB.
"""
function lcdb_entry(object::AbstractString)
    page, _ = _alcdef_objectpage(object)
    i = findfirst("<h3>Summary</h3>", page)
    i === nothing && return Dict{String,String}()
    tab = match(r"<table.*?</table>"s, page[last(i):end])
    tab === nothing && return Dict{String,String}()
    rows = [m.captures[1] for m in eachmatch(r"<tr[^>]*>(.*?)</tr>"s, tab.match)]
    length(rows) < 2 && return Dict{String,String}()
    hdr = [_html_text(m.captures[1]) for m in eachmatch(r"<th[^>]*>(.*?)</th>"s, rows[1])]
    val = [_html_text(m.captures[1]) for m in eachmatch(r"<td[^>]*>(.*?)</td>"s, rows[2])]
    out = Dict{String,String}()
    for (h, v) in zip(hdr, val)
        haskey(out, h) && (h = h * "2")              # "AF" appears twice (albedo flag, amplitude flag)
        out[h] = v
    end
    return out
end

const LCDB_ZIP_URL = "https://minplanobs.org/MPInfo/datazips/LCLIST_PUB_CURRENT.zip"

"""
    fetch_lcdb_member(name="lc_summary_pub.txt"; url=LCDB_ZIP_URL, dir=datadir(), force=false) -> path

Download one file of the LCDB public release zip (≈ 40 MB in total) by HTTP range requests: only the zip directory and
the requested member are transferred (`lc_summary_pub.txt`: ≈ 1.3 MB compressed, 8 MB uncompressed); if the whole zip
is in `dir` ([`fetch_lcdb_archive`](@ref)), the member is extracted locally. Cached in `dir`.
Other members: `lc_details_pub.txt` (≈ 4.8 MB compressed), `lc_spinaxis_pub.txt`, `lc_binary_pub.txt`, … .
"""
function fetch_lcdb_member(name::AbstractString="lc_summary_pub.txt"; url::AbstractString=LCDB_ZIP_URL,
                           dir::AbstractString=datadir(), force::Bool=false)
    path = joinpath(dir, name)
    isfile(path) && !force && return path
    mkpath(dir)
    localzip = joinpath(dir, "LCLIST_PUB_CURRENT.zip")
    write(path, _zip_member(isfile(localzip) ? localzip : url, name))
    return path
end

# zip access for a URL (HTTP range requests) or a local file: entries of the central directory and single members
function _zipsource(src::AbstractString)
    if isfile(src)
        sz = filesize(src)
        rd = (a, b) -> open(io -> (seek(io, a); read(io, b - a + 1)), src)
        return sz, rd
    end
    _, resp = _http(src; method="HEAD")
    hs = Dict(lowercase(k) => v for (k, v) in resp.headers)
    return parse(Int, hs["content-length"]), (a, b) -> first(_http(src; headers=["Range" => "bytes=$a-$b"]))
end
_u16(v, i) = UInt16(v[i]) | UInt16(v[i+1]) << 8
_u32(v, i) = UInt32(v[i]) | UInt32(v[i+1]) << 8 | UInt32(v[i+2]) << 16 | UInt32(v[i+3]) << 24

function _zip_entries(sz, rd)
    tail = rd(max(sz - 65557, 0), sz - 1)
    e = findlast(i -> _u32(tail, i) == 0x06054b50, 1:length(tail)-3)
    e === nothing && throw(ErrorException("zip: end of central directory not found"))
    cdsize = Int(_u32(tail, e + 12)); cdoff = Int(_u32(tail, e + 16))
    cd = rd(cdoff, cdoff + cdsize - 1)
    out = NamedTuple{(:name, :method, :csize, :usize, :offset),Tuple{String,Int,Int,Int,Int}}[]
    i = 1
    while i + 46 <= length(cd) && _u32(cd, i) == 0x02014b50
        nlen = _u16(cd, i + 28); xlen = _u16(cd, i + 30); clen = _u16(cd, i + 32)
        push!(out, (name=String(cd[i+46:i+45+nlen]), method=Int(_u16(cd, i + 10)), csize=Int(_u32(cd, i + 20)),
                    usize=Int(_u32(cd, i + 24)), offset=Int(_u32(cd, i + 42))))
        i += 46 + nlen + xlen + clen
    end
    return out
end
function _zip_read(rd, e)
    lh = rd(e.offset, e.offset + 29)
    start = e.offset + 30 + _u16(lh, 27) + _u16(lh, 29)
    data = rd(start, start + e.csize - 1)
    out = e.method == 0 ? data : e.method == 8 ? transcode(CodecZlib.DeflateDecompressor, data) :
          throw(ErrorException("zip: compression method $(e.method) not supported"))
    length(out) == e.usize || throw(ErrorException("zip: size mismatch for $(e.name)"))
    return out
end
function _zip_member(src, name)
    sz, rd = _zipsource(src)
    i = findfirst(e -> e.name == name, _zip_entries(sz, rd))
    i === nothing && throw(ErrorException("zip: member $name not found in $src"))
    return _zip_read(rd, _zip_entries(sz, rd)[i])
end

function _download(url, path; force=false)
    isfile(path) && !force && return path
    mkpath(dirname(path))
    tmp = path * ".part"
    Downloads.download(url, tmp; headers=["User-Agent" => USER_AGENT])
    mv(tmp, path; force=true)
    return path
end

"""
    fetch_lcdb_archive(; dir=datadir(), force=false) -> path

Download the complete LCDB public release (`LCLIST_PUB_CURRENT.zip`, ≈ 40 MB) to `dir`. [`fetch_lcdb_member`](@ref)
then reads its members locally.
"""
fetch_lcdb_archive(; dir::AbstractString=datadir(), force::Bool=false) =
    _download(LCDB_ZIP_URL, joinpath(dir, "LCLIST_PUB_CURRENT.zip"); force=force)

const ALCDEF_ARCHIVE_URL = "https://www.alcdef.org/docs/ALCDEF_ALL.zip"

"""
    fetch_alcdef_archive(; dir=datadir(), force=false) -> path

Download the complete ALCDEF archive (`ALCDEF_ALL.zip`, ≈ 135 MB, one file per object, updated about every three months)
to `dir`. Read objects from it with [`read_alcdef_archive`](@ref).
"""
fetch_alcdef_archive(; dir::AbstractString=datadir(), force::Bool=false) =
    _download(ALCDEF_ARCHIVE_URL, joinpath(dir, "ALCDEF_ALL.zip"); force=force)

"""
    alcdef_archive_objects(archive=joinpath(datadir(), "ALCDEF_ALL.zip")) -> Vector{NamedTuple}

The objects in the local ALCDEF archive: `(number, name, file, bytes)` from the member names
`ALCDEF_<number>_<name>.txt` (number 0 for unnumbered objects, whose name is the designation).
"""
function alcdef_archive_objects(archive::AbstractString=joinpath(datadir(), "ALCDEF_ALL.zip"))
    isfile(archive) || throw(ArgumentError("$archive not found; download it with fetch_alcdef_archive()"))
    sz, rd = _zipsource(archive)
    out = NamedTuple{(:number, :name, :file, :bytes),Tuple{Int,String,String,Int}}[]
    for e in _zip_entries(sz, rd)
        m = match(r"^ALCDEF_(\d+)_(.+)\.txt$", e.name)
        m === nothing && continue
        push!(out, (number=parse(Int, m.captures[1]), name=replace(m.captures[2], "_" => " "), file=e.name, bytes=e.usize))
    end
    return out
end

"""
    read_alcdef_archive(object; archive=joinpath(datadir(), "ALCDEF_ALL.zip")) -> Vector{ALCDEFBlock}

The ALCDEF light curves of one object from the local archive ([`fetch_alcdef_archive`](@ref)), by number (Integer or
digits) or name / designation (case-insensitive; spaces or underscores). Empty if the object is not in the archive.
"""
function read_alcdef_archive(object; archive::AbstractString=joinpath(datadir(), "ALCDEF_ALL.zip"))
    objs = alcdef_archive_objects(archive)
    key = strip(string(object))
    i = if all(isdigit, key) && key != "0"
        findfirst(o -> o.number == parse(Int, key), objs)
    else
        k = lowercase(replace(key, "_" => " "))
        findfirst(o -> lowercase(o.name) == k, objs)
    end
    i === nothing && return ALCDEFBlock[]
    sz, rd = _zipsource(archive)
    e = only(filter(x -> x.name == objs[i].file, _zip_entries(sz, rd)))
    return read_alcdef(IOBuffer(_zip_read(rd, e)))
end

"""
    read_lcdb_summary(path=fetch_lcdb_member()) -> NamedTuple

Read the LCDB summary table (`lc_summary_pub.txt`, fixed-width columns located from its header line). Returns the
columns as `String` vectors keyed by the header names (e.g. `tab.columns["PERIOD"]`) plus typed vectors `number`,
`name`, `desig`, `diameter`, `H`, `albedo`, `period` [h], `amin`, `amax`, `U` (NaN / "" where empty).
"""
function read_lcdb_summary(path::AbstractString=fetch_lcdb_member())
    L = readlines(path; keep=false)
    ih = findfirst(l -> startswith(l, "NUMBER"), L)
    ih === nothing && throw(ArgumentError("$path: header line not found"))
    h = L[ih]
    toks = [(m.offset, m.match) for m in eachmatch(r"\S+", h)]
    names = String[]; seen = Dict{String,Int}()
    for (_, t) in toks
        n = get(seen, t, 0); seen[t] = n + 1
        push!(names, n == 0 ? t : "$(t)_$(n+1)")      # repeated flag columns S, F: S, S_2, …
    end
    starts = first.(toks)
    data = filter(l -> !isempty(strip(l)) && !startswith(l, "---"), L[ih+1:end])
    cols = Dict(n => String[] for n in names)
    for l in data
        for (j, n) in enumerate(names)
            a = starts[j]; b = j < length(starts) ? starts[j+1] - 1 : length(l)
            push!(cols[n], a > length(l) ? "" : strip(l[a:min(b, length(l))]))
        end
    end
    f(n) = [something(tryparse(Float64, x), NaN) for x in cols[n]]
    return (columns=cols, number=[something(tryparse(Int, x), 0) for x in cols["NUMBER"]], name=cols["NAME"],
            desig=cols["DESIG"], diameter=f("DIA."), H=f("H"), albedo=f("ALBEDO"), period=f("PERIOD"), amin=f("AMIN"),
            amax=f("AMAX"), U=cols["U"])
end

"""
    lcdb_lookup(tab, key) -> Int or nothing

Row index of an object in a table from [`read_lcdb_summary`](@ref), by number (Integer) or name/designation (String).
"""
lcdb_lookup(tab, key::Integer) = findfirst(==(key), tab.number)
lcdb_lookup(tab, key::AbstractString) = findfirst(i -> lowercase(tab.name[i]) == lowercase(key) ||
                                                       lowercase(tab.desig[i]) == lowercase(key), eachindex(tab.name))
