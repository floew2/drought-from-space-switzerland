"""The sensor response and pixel purity. METHOD §2.3.

Purity is the share of a 10 m cell's sensor RESPONSE that comes from inside a
category's fields. Fields are rasterised at 2.5 m, convolved with the modelled
response, and read per cell. The response is the detector box convolved with
Gaussian optics; the box is applied as a full-cell average.
"""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from PIL import Image, ImageDraw

from .settings import PSF_MODEL, PSF_MTF, RES, SUBPIXEL

__all__ = ["measured_response_fwhm", "psf_caption", "psf_params", "psf_sigma_m",
           "purity_psf"]


def psf_params(res: float | None = None, mtf: float | None = None,
               model: str | None = None) -> dict:
    """The modelled sensor response, derived from the system MTF at Nyquist.

    For a Gaussian, MTF(f) = exp(-2 pi^2 sigma^2 f^2) with f_Nyquist = 1/(2p).

    box_gauss  The detector footprint is a box one pixel wide, with its own MTF
               sinc(1/2) = 0.637 at Nyquist. Only the remainder of the system MTF
               belongs to the optics: 0.25 / 0.637 = 0.393, so sigma_optics =
               4.35 m. The combined response (box convolved with that Gaussian)
               has FWHM 12.7 m; the Gaussian alone would be 10.2 m.
    gauss      One Gaussian fitted to the whole system MTF: sigma 5.30 m, FWHM
               12.5 m. The earlier implementation, which built the mask behind
               the published area numbers; kept for comparison. As sampled
               (a 2 x 2 sub-cell mean, not a point) its effective response is
               ~12.95 m — see `measured_response_fwhm`.

    Args:
        res: Pixel size in metres (default: config).
        mtf: System MTF at Nyquist (default: config).
        model: "box_gauss" or "gauss" (default: config).

    Returns:
        dict with model, sigma_m (the Gaussian term), fwhm_m (of the whole
        response), mtf_system, mtf_detector, mtf_optics.

    Raises:
        ValueError: For an unknown model.
    """
    from scipy.optimize import brentq
    from scipy.special import erf

    p_ = RES if res is None else res
    m_ = PSF_MTF if mtf is None else mtf
    model = PSF_MODEL if model is None else model
    f_n = 1.0 / (2.0 * p_)
    if model == "box_gauss":
        mtf_det = float(np.sinc(p_ * f_n))           # numpy sinc is sin(pi x)/(pi x)
    elif model == "gauss":
        mtf_det = 1.0                                # detector folded into the fit
    else:
        raise ValueError(f"unknown psf_model {model!r}; use 'box_gauss' or 'gauss'")
    mtf_opt = m_ / mtf_det
    sigma = float(np.sqrt(-np.log(mtf_opt) / (2.0 * np.pi**2 * f_n**2)))
    if model == "box_gauss":
        def prof(x):
            s2 = sigma * np.sqrt(2.0)
            return 0.5 * (erf((x + p_ / 2) / s2) - erf((x - p_ / 2) / s2))
        half = brentq(lambda x: prof(x) / prof(0.0) - 0.5, 0.0, 10 * p_)
        fwhm = 2.0 * half
    else:
        fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * sigma
    return {"model": model, "sigma_m": sigma, "fwhm_m": float(fwhm),
            "mtf_system": float(m_), "mtf_detector": mtf_det,
            "mtf_optics": float(mtf_opt)}


def psf_sigma_m(res: float | None = None, mtf: float | None = None,
                model: str | None = None) -> float:
    """Sigma of the Gaussian term of the modelled response, in metres.

    Under the default box_gauss model this is the OPTICAL sigma (4.35 m), not
    the width of the whole response; `psf_params` gives that.

    Args:
        res: Pixel size in metres (default: config).
        mtf: System MTF at Nyquist (default: config).
        model: "box_gauss" or "gauss" (default: config).

    Returns:
        Sigma in metres.
    """
    return psf_params(res, mtf, model)["sigma_m"]


def measured_response_fwhm(model: str | None = None, sigma_m: float | None = None) -> float:
    """FWHM of the response that `purity_psf` actually applies, measured.

    `psf_params` states what the model should be; this measures what the code
    does, by passing a line source (one 2.5 m sub-cell wide) through
    `purity_psf` at four sub-cell offsets and reading the profile at 2.5 m
    spacing. The line source itself adds its 2.5 m width, so the expected
    values are box(2.5) (x) response:

        box_gauss   box(10) (x) G(4.35)            12.67 m -> ~12.8 measured
        gauss       2x2 sub-cell mean (x) G(5.30)  ~12.95 m -> ~13.1 measured
                    (the four-point "centre" estimate spans 5 m, so the old
                    implementation was wider than its nominal 12.48 m)
        gauss, sigma 4.35 — swapping sigma without averaging the cell:
                                                   ~10.8 m -> ~11.0 measured

    Args:
        model: "box_gauss" or "gauss" (default: config).
        sigma_m: Override for the Gaussian sigma.

    Returns:
        Measured FWHM in metres.
    """
    from shapely.geometry import box as _box

    n, res, xs, vs = 21, float(RES), [], []
    for k in range(SUBPIXEL):
        x0 = 10 * res + k * res / SUBPIXEL
        p = purity_psf([_box(x0, 0, x0 + res / SUBPIXEL, n * res)],
                       (0.0, n * res, res), (n, n), model=model, sigma_m=sigma_m)
        for c in range(n):
            xs.append((c + 0.5) * res - (x0 + res / SUBPIXEL / 2))
            vs.append(float(p[n // 2, c]))
    o = np.argsort(xs)
    x, v = np.asarray(xs)[o], np.asarray(vs)[o]
    v = v / v.max()
    i0, i1 = np.where(v >= 0.5)[0][[0, -1]]
    left = np.interp(0.5, [v[i0 - 1], v[i0]], [x[i0 - 1], x[i0]])
    right = np.interp(0.5, [v[i1 + 1], v[i1]], [x[i1 + 1], x[i1]])
    return float(right - left)


def psf_caption(model: str | None = None) -> str:
    """One-line description of the response model, for figure captions.

    Args:
        model: "box_gauss" or "gauss" (default: config).

    Returns:
        Caption text.
    """
    p = psf_params(model=model)
    if p["model"] == "box_gauss":
        return (f"{RES:g} m detector box \u2297 Gaussian optics (\u03c3 {p['sigma_m']:.2f} m), "
                f"FWHM {p['fwhm_m']:.1f} m, from system MTF {p['mtf_system']:g} at Nyquist")
    return (f"Gaussian \u03c3 {p['sigma_m']:.2f} m, FWHM {p['fwhm_m']:.1f} m, "
            f"from system MTF {p['mtf_system']:g} at Nyquist")


def _burn(geoms, grid, shape) -> np.ndarray:
    """Rasterise polygons: 1 where a cell centre falls inside a polygon.

    Uses GDAL's rasteriser via rasterio, whose centre-in-polygon rule does not
    depend on where the window starts — so a tiled build gives exactly the
    single-pass result. PIL's polygon fill, used earlier to sidestep a
    rasterio/affine version conflict, flips edge pixels depending on the window
    offset (about 50 edge pixels in a 172 x 172 window), which made tiled purity
    differ from untiled purity by up to 0.15. It remains only as a fallback.

    Args:
        geoms: Iterable of shapely (Multi)Polygons.
        grid: (minx, maxy, res) of the raster.
        shape: (height, width) in cells.

    Returns:
        uint8 array of 0/1.
    """
    geoms = [g for g in geoms if g is not None and not g.is_empty]
    if not geoms:
        return np.zeros(shape, np.uint8)
    minx, maxy, res = grid
    try:
        from rasterio import features
        from rasterio.transform import from_origin
    except ImportError:                                   # pragma: no cover
        return _burn_pil(geoms, grid, shape)
    return features.rasterize(((g, 1) for g in geoms), out_shape=shape,
                              transform=from_origin(minx, maxy, res, res),
                              fill=0, dtype="uint8")


def _burn_pil(geoms, grid, shape) -> np.ndarray:
    """Fallback rasteriser (PIL). Edge pixels depend on the window offset."""
    minx, maxy, res = grid
    img = Image.new("1", (shape[1], shape[0]), 0)
    draw = ImageDraw.Draw(img)
    for geom in geoms:
        if geom is None or geom.is_empty:
            continue
        parts = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for poly in parts:
            xs, ys = poly.exterior.coords.xy
            px = [(x - minx) / res for x in xs]
            py = [(maxy - y) / res for y in ys]
            if len(px) > 2:
                draw.polygon(list(zip(px, py)), fill=1)
            for ring in poly.interiors:
                xs, ys = ring.coords.xy
                px = [(x - minx) / res for x in xs]
                py = [(maxy - y) / res for y in ys]
                if len(px) > 2:
                    draw.polygon(list(zip(px, py)), fill=0)
    return np.array(img, dtype=np.uint8)


def purity_psf(geoms: Iterable, grid: tuple[float, float, float], shape: tuple[int, int],
               sub_n: int | None = None, sigma_m: float | None = None,
               model: str | None = None, tile: int = 512) -> np.ndarray:
    """PSF-weighted purity: the share of a cell's RESPONSE from inside a parcel.

    Rasterise at res/sub_n (2.5 m), convolve with the modelled response, and
    read the value for each 10 m cell. What contaminates the reflectance is the
    response, not the area, so this is the quantity a plain area fraction only
    approximates.

    Under box_gauss the Gaussian (optics) is applied on the fine grid and the
    detector box is applied by averaging the sub_n x sub_n sub-cells of each
    cell — which is exactly the box convolution evaluated at the cell centre.
    Under gauss the single Gaussian is sampled at the cell centre.

    Computed in tiles with a halo wider than the response, so memory stays
    bounded for a canton-sized grid and the result is identical to one pass.

    Args:
        geoms: Iterable of shapely (Multi)Polygons in the grid CRS.
        grid: (minx, maxy, res) of the output grid.
        shape: (height, width) of the output grid in cells.
        sub_n: Sub-cells per side (default: config, 4 -> 2.5 m).
        sigma_m: Override for the Gaussian sigma in metres.
        model: "box_gauss" or "gauss" (default: config).
        tile: Output tile size in cells.

    Returns:
        float32 array of purity in [0, 1], shape `shape`.
    """
    from scipy.ndimage import gaussian_filter
    from shapely.geometry import box as _box
    from shapely.strtree import STRtree

    sub_n = SUBPIXEL if sub_n is None else sub_n
    model = PSF_MODEL if model is None else model
    minx, maxy, res = grid
    height, width = shape
    fine_res = res / sub_n
    sigma = psf_sigma_m(model=model) if sigma_m is None else sigma_m
    sigma_px = sigma / fine_res
    # gaussian_filter truncates at 4 sigma; the halo must cover that plus the
    # half-cell of the detector box, so tile edges never see a missing parcel.
    halo = int(np.ceil((4.0 * sigma_px + 1) / sub_n)) + 1

    geoms = [g for g in geoms if g is not None and not g.is_empty]
    out = np.zeros((height, width), np.float32)
    if not geoms:
        return out
    tree = STRtree(geoms)

    for r0 in range(0, height, tile):
        r1 = min(r0 + tile, height)
        for c0 in range(0, width, tile):
            c1 = min(c0 + tile, width)
            R0, R1, C0, C1 = r0 - halo, r1 + halo, c0 - halo, c1 + halo
            x0, x1 = minx + C0 * res, minx + C1 * res
            y1, y0 = maxy - R0 * res, maxy - R1 * res
            hit = tree.query(_box(x0, y0, x1, y1))
            if not len(hit):
                continue
            nr, nc = R1 - R0, C1 - C0
            a = _burn([geoms[i] for i in hit], (x0, y1, fine_res),
                      (nr * sub_n, nc * sub_n)).astype(np.float32)
            a = gaussian_filter(a, sigma_px, mode="constant", cval=0.0)
            if model == "box_gauss":
                cell = a.reshape(nr, sub_n, nc, sub_n).mean(axis=(1, 3))
            elif sub_n % 2 == 0:             # the centre falls between sub-cells
                c = sub_n // 2
                cell = 0.25 * (a[c - 1::sub_n, c - 1::sub_n] + a[c - 1::sub_n, c::sub_n]
                               + a[c::sub_n, c - 1::sub_n] + a[c::sub_n, c::sub_n])
            else:
                cell = a[sub_n // 2::sub_n, sub_n // 2::sub_n]
            out[r0:r1, c0:c1] = cell[halo:halo + (r1 - r0), halo:halo + (c1 - c0)]
    # float32 convolution leaves a fully interior cell at 0.99999994, which
    # would make "purity >= 1" select nothing. Six decimals is far below any
    # threshold anyone sets.
    return np.clip(np.round(out, 6), 0.0, 1.0).astype(np.float32)
