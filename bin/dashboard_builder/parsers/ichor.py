"""
ichor.py - dashboard parser for ichorCNA output.

The copy-number tab shows the ichorCNA genome-wide figure and the fitted
parameters. It deliberately does not show a segment call table: large-scale
copy number is read off the plot, and the per-bin segment file is not a
clinical reporting artefact.

Inputs, as laid down by build_report_bundle.sh:

    <effective_dir>/cnv/<sample>.ichor_all_sols.pdf
    <effective_dir>/cnv/<sample>.ichor_params.txt

Rendering strategy
------------------
A PDF referenced from an <iframe> or <embed> renders inconsistently when the
report is opened over file://, and Chrome blocks navigation to a base64 data:
URI containing a PDF. So the pages are rasterised once at build time with
pdftoppm and inlined as base64 PNGs, which display everywhere. The original
PDF is still copied into the bundle and linked, so the vector version remains
one click away for anyone who wants to zoom.

If pdftoppm is unavailable the parser degrades to an <object> embed plus the
download link rather than failing the build.

Standard library only.
"""

import base64
import glob
import os
import re
import shutil
import subprocess
import tempfile


# Rasterisation settings. 110 dpi keeps a genome-wide ichorCNA panel legible
# at full width while holding a typical page under ~250 kB of base64.
RENDER_DPI = 110
RENDER_MAX_WIDTH = 1600

# Guard against inlining an unbounded number of solutions. The all-solutions
# PDF carries one page per (normal fraction, ploidy) combination tried.
MAX_INLINE_PAGES = 16


def _find_pdf(effective_dir, sample):
    """Locate the ichorCNA figure PDF, preferring the all-solutions panel."""
    direct = os.path.join(
        effective_dir, "cnv", "%s.ichor_all_sols.pdf" % sample
    )
    if os.path.isfile(direct):
        return direct

    candidates = []
    for root, _dirs, files in os.walk(effective_dir):
        for name in files:
            if not name.lower().endswith(".pdf"):
                continue
            if sample not in name and sample not in root:
                continue
            lowered = name.lower()
            # Rank: all-solutions panel first, then any genome-wide figure.
            if "all_sols" in lowered or "all_sol" in lowered:
                rank = 0
            elif "genomewide" in lowered.replace("_", ""):
                rank = 1
            elif "ichor" in lowered or "ichor" in root.lower():
                rank = 2
            else:
                continue
            candidates.append((rank, name, os.path.join(root, name)))

    if not candidates:
        return None
    candidates.sort()
    return candidates[0][2]


def _find_params(effective_dir, sample):
    """Locate the ichorCNA params.txt for this sample."""
    direct = os.path.join(
        effective_dir, "cnv", "%s.ichor_params.txt" % sample
    )
    if os.path.isfile(direct):
        return direct

    for root, _dirs, files in os.walk(effective_dir):
        for name in files:
            if not name.endswith("params.txt"):
                continue
            if sample in name or sample in root:
                return os.path.join(root, name)
    return None


def _parse_params(path):
    """Parse ichorCNA params.txt into an ordered list of key/value pairs.

    The file is a small tab-delimited block; ichorCNA has emitted it both as a
    two-line header/value table and as one key-value pair per line across
    versions. Both shapes are handled by inspecting the first two lines.
    """
    if not path or not os.path.isfile(path):
        return []

    with open(path) as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    if not lines:
        return []

    first = lines[0].split("\t")
    second = lines[1].split("\t") if len(lines) > 1 else []

    # Header/value table: two lines of equal, greater-than-two width.
    if len(first) > 2 and len(second) == len(first):
        return [
            {"key": k.strip(), "value": v.strip()}
            for k, v in zip(first, second)
            if k.strip()
        ]

    # Key-value per line.
    pairs = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) >= 2:
            pairs.append({"key": parts[0].strip(), "value": parts[1].strip()})
        elif ":" in line:
            key, value = line.split(":", 1)
            pairs.append({"key": key.strip(), "value": value.strip()})
    return pairs


def _headline(params):
    """Pull the two numbers a reader looks for first, if they are present.

    Matching is on the parameter name only, and returns whatever the fit
    produced. No expected or reference value is encoded.
    """
    tumour = None
    ploidy = None
    for pair in params:
        key = pair["key"].lower().replace("_", " ").strip()
        if tumour is None and (
            "tumor fraction" in key or "tumour fraction" in key
        ):
            tumour = pair["value"]
        if ploidy is None and key == "ploidy":
            ploidy = pair["value"]
    return {"tumour_fraction": tumour, "ploidy": ploidy}


def _page_sort_key(path):
    """Sort rasterised pages numerically (page-2 before page-10)."""
    match = re.search(r"-(\d+)\.png$", path)
    return int(match.group(1)) if match else 0


def _render_pages(pdf_path, max_pages=MAX_INLINE_PAGES):
    """Rasterise a PDF to a list of base64 PNG data URIs.

    Returns an empty list when pdftoppm is not installed or the conversion
    fails, which the template treats as "fall back to the embedded object".
    """
    if not shutil.which("pdftoppm"):
        return []

    workdir = tempfile.mkdtemp(prefix="ichor_png_")
    try:
        command = [
            "pdftoppm",
            "-png",
            "-r",
            str(RENDER_DPI),
            "-scale-to-x",
            str(RENDER_MAX_WIDTH),
            "-scale-to-y",
            "-1",
            "-f",
            "1",
            "-l",
            str(max_pages),
            pdf_path,
            os.path.join(workdir, "page"),
        ]
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        if result.returncode != 0:
            return []

        pages = sorted(
            glob.glob(os.path.join(workdir, "page-*.png")), key=_page_sort_key
        )
        images = []
        for page in pages[:max_pages]:
            with open(page, "rb") as handle:
                encoded = base64.b64encode(handle.read()).decode("ascii")
            images.append("data:image/png;base64," + encoded)
        return images
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def parse(effective_dir, sample, max_pages=MAX_INLINE_PAGES):
    """Collect the ichorCNA figure and fit parameters for one sample."""
    pdf_path = _find_pdf(effective_dir, sample)
    params_path = _find_params(effective_dir, sample)

    if not pdf_path and not params_path:
        return {
            "found": False,
            "reason": "no ichorCNA PDF or params.txt for this sample",
            "searched": str(effective_dir),
            "images": [],
            "params": [],
            "headline": {"tumour_fraction": None, "ploidy": None},
            "n_pages": 0,
            "rasterised": False,
            "pdf_href": None,
        }

    params = _parse_params(params_path)
    images = _render_pages(pdf_path, max_pages) if pdf_path else []

    # Relative href from the sample report HTML, which is written into the
    # sample directory itself.
    relative_pdf = None
    if pdf_path:
        try:
            relative_pdf = os.path.relpath(pdf_path, effective_dir)
        except ValueError:
            relative_pdf = None

    return {
        "found": True,
        "searched": str(effective_dir),
        "pdf_path": pdf_path,
        "pdf_href": relative_pdf,
        "pdf_filename": os.path.basename(pdf_path) if pdf_path else None,
        "params_path": params_path,
        "params": params,
        "headline": _headline(params),
        "images": images,
        "n_pages": len(images),
        "rasterised": bool(images),
    }
