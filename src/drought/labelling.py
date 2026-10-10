"""Interactive labelling widget.

Pure HTML/JS embedded in the notebook — no ipympl, no plotly, nothing to
install, and no kernel round-trip, so it cannot crash the cell.

Reference data is per SITE: a detector validated on one site has not been
validated on another, and pooling the labels would hide that.
"""

from __future__ import annotations

import json

import numpy as np

from .config import CFG, labels_path
from .detect import regime_of


def label_widget(df, site, pairs=None, n=20, seed=0, ev=None,
                 codes=None, height=300, existing=None,
                 skip_labelled=True):
    """Click-to-label tool for grassland cut events.

    Pure HTML/JS embedded in the notebook — no ipympl, no plotly, nothing to
    install, and no kernel round-trip, so it cannot crash the cell.

    Controls
        click on the plot      add a label at the nearest observation
        click a blue marker    remove that label
        A                      accept the detector's cuts as labels
        H                      show / hide the detector's proposals
        C                      clear labels for this parcel-year
        <- / ->                previous / next parcel-year

    The detector's proposals are faint red dashed lines. Press H to hide them
    and label blind, then H again to compare — a label set that agrees with the
    detector because it was copied from it measures nothing. Blind-then-reveal
    is the honest workflow; A is the fast one.

    Place each label at the TROUGH: the first point that is clearly down, not
    the midpoint and not the last high point. Clicks snap to the nearest
    observation, so a label always sits on real data.

    Each panel is headed with its LNF code and regime, because the two behave
    completely differently: 613 takes 4-6 cuts from early May, while 611/612
    take about two and may not legally be cut before 15 June — shown as a
    shaded zone so a label placed there is visibly questionable.

    Previous work is picked up automatically from `labels_path`
    (default cut_labels.json) when `existing` is not given, and those
    parcel-years are excluded from the new sample. Pass `existing={}` to
    deliberately start fresh, or `existing=` a dict to continue: those parcel-years are
    skipped in the sample and their labels are carried into the output, so the
    box always holds the COMPLETE set rather than just this session's work.

    Output accumulates in the text box at the bottom in exactly the format
    `score()` and `tune()` expect. Copy it into a cell as `labels = {...}`.
    """
    from IPython.display import HTML, display

    # Default to every mown code, not just 613: the regime matters for
    # labelling — an extensive meadow (611/612) legally cannot be cut before
    # 15 June and takes ~2 cuts, while 613 takes 4-6 from early May. Labelling
    # one as if it were the other produces bad ground truth.
    if codes is None:
        codes = list(CFG.landuse.intensive) + list(CFG.landuse.extensive)
    d = df[df.lnf_code.isin(codes)] if "lnf_code" in df else df

    # Pick up previous work automatically. Requiring the caller to remember
    # `existing=` is how a labelling session silently starts from scratch and
    # re-proposes parcel-years that were already done.
    #   existing=None  -> load labels_path if it exists
    #   existing={...} -> use exactly that
    #   existing={}    -> start fresh on purpose
    if existing is None:
        lp = labels_path(site)
        if lp.exists():
            try:
                with open(lp) as f:
                    existing = json.load(f)
                nc = sum(len(v) for v in existing.values())
                print(f"resuming from {lp}: {len(existing)} "
                      f"parcel-years, {nc} cuts already labelled")
            except Exception as exc:
                print(f"  [!] could not read {lp} ({exc}) — "
                      "starting with no prior labels")
                existing = {}
        else:
            existing = {}
    existing = dict(existing or {})
    rng = np.random.default_rng(seed)
    if pairs is None:
        combos = d.groupby(["parcel_id", "year"]).size()
        pairs = combos[combos >= 12].index.tolist()

    # `n` now applies whether or not `pairs` was supplied. It previously only
    # applied to the random sample, so passing an explicit list built a panel
    # for EVERY parcel-year — megabytes of embedded JSON, which Jupyter
    # truncates, cutting the script mid-array so nothing renders at all.
    pairs = [tuple(x) for x in pairs]
    if skip_labelled and existing:
        done = set(existing)
        pairs = [c for c in pairs if f"{c[0]}|{int(c[1])}" not in done]
        if not pairs:
            print("every candidate parcel-year is already labelled — widen "
                  "`codes`/`pairs`, or pass `existing={}` to relabel")
            return
    print(f"{len(pairs)} unlabelled parcel-years available")

    if n and len(pairs) > n:
        idx = rng.choice(len(pairs), size=n, replace=False)
        pairs = [pairs[i] for i in sorted(idx)]

    panels = []
    for pid, yr in pairs:
        g = d[(d.parcel_id == pid) & (d.year == yr)].sort_values("doy")
        if g.empty:
            continue
        haze = g["haze"].tolist() if "haze" in g else [False] * len(g)
        cuts = []
        if ev is not None:
            r = ev[(ev.parcel_id == pid) & (ev.year == yr)]
            if len(r) and isinstance(r.iloc[0].get("cut_dates"), list):
                cuts = [float(c) for c in r.iloc[0]["cut_dates"]]
        code = int(g["lnf_code"].iloc[0]) if "lnf_code" in g else None
        reg = regime_of(code) if code is not None else "unknown"
        panels.append({
            "key": f"{pid}|{int(yr)}",
            "pid": str(pid)[-8:],
            "pid_full": str(pid),
            "year": int(yr),
            "code": code,
            "regime": reg,
            "expect": {"intensive": "4-6 cuts, first early May, no legal limit",
                       "extensive": "~2 cuts, none before 15 June (DOY 166)"}
                      .get(reg, ""),
            "legal": 166 if reg == "extensive" else None,
            "doy": [int(v) for v in g["doy"]],
            "ndmi": [round(float(v), 4) for v in g["NDMI"]],
            "ndvi": [round(float(v), 4) for v in g["NDVI"]],
            "haze": [bool(h) for h in haze],
            "cuts": cuts,
            "have": [int(x) for x in existing.get(f"{pid}|{int(yr)}", [])],
        })

    # The header truncates the id for space. Make sure that truncation cannot
    # be ambiguous, since the label key uses the FULL id and a reader matching
    # screen to file would otherwise be misled.
    short = {}
    for p in panels:
        short.setdefault(p["pid"], set()).add(p["pid_full"])
    clash = {k: v for k, v in short.items() if len(v) > 1}
    if clash:
        print(f"  [!] truncated ids are ambiguous for {list(clash)[:3]} — "
              "read the full id shown under the header")

    if not panels:
        print("no parcel-years matched — check that `pairs` is not empty and "
              "that its parcel ids and years exist in the data")
        return
    thin = [p["key"] for p in panels if len(p["doy"]) < 5]
    if thin:
        print(f"  [!] {len(thin)} panel(s) have fewer than 5 observations and "
              f"will look empty: {thin[:5]}")

    # Labels already on disk are carried into the widget, so the output box
    # always holds the COMPLETE set — never just this session's additions.
    prior = {k: [int(x) for x in v] for k, v in existing.items()}
    payload = json.dumps(panels)
    if len(payload) > 1_500_000:
        print(f"  [!] {len(payload)/1e6:.1f} MB of data for {len(panels)} "
              "panels. Jupyter truncates large output, which breaks the widget "
              "— lower `n`.")
    html = (_LABEL_HTML.replace("__DATA__", payload)
                       .replace("__PRIOR__", json.dumps(prior))
                       .replace("__HEIGHT__", str(height)))
    display(HTML(html))
    print(f"{len(panels)} parcel-years loaded. Label them, then copy the dict "
          "from the box at the bottom into a cell as:  labels = {...}")


_LABEL_HTML = r"""
<div id="lblwrap" style="font-family:system-ui,sans-serif;max-width:1100px">
  <div style="display:flex;align-items:center;gap:10px;margin-bottom:6px">
    <button id="prev">&larr; prev</button>
    <button id="next">next &rarr;</button>
    <span id="pos" style="font-weight:600;line-height:1.25"></span>
    <span id="cnt" style="color:#666;font-size:12px"></span>
    <span style="flex:1"></span>
    <button id="hide" title="show / hide detector proposals (H)">hide detector</button>
    <button id="accept" title="accept detector cuts (A)">accept detector</button>
    <button id="clear" title="clear this panel (C)">clear</button>
  </div>
  <svg id="plot" width="1080" height="__HEIGHT__"
       style="border:1px solid #ccc;background:#fff;cursor:crosshair"></svg>
  <div style="font-size:12px;color:#555;margin:4px 0">
    click = add label at nearest observation &nbsp;|&nbsp; click a blue marker =
    remove &nbsp;|&nbsp; blue = your labels, red dashed = detector
  </div>
  <div id="err" style="display:none;color:#b00;font-size:12px;
       font-family:ui-monospace,monospace;padding:6px;background:#fee;
       border:1px solid #fbb;margin:4px 0"></div>
  <textarea id="out" style="width:100%;height:130px;font-family:ui-monospace,
    monospace;font-size:12px"></textarea>
</div>
<script>
(function(){
  // Without this a single JS exception leaves an empty SVG and no explanation.
  function fail(m){
    const e = document.getElementById('err');
    e.style.display = 'block';
    e.textContent = 'widget error: ' + m;
  }
  window.addEventListener('error', ev => fail(ev.message));
  try {
  const DATA = __DATA__;
  const H = __HEIGHT__, W = 1080;
  const M = {l:48, r:14, t:26, b:30};
  let X0 = 95, X1 = 300, Y0 = 0, Y1 = 1;
  (function fit(){
    let lo = 1e9, hi = -1e9, ylo = 1e9, yhi = -1e9;
    DATA.forEach(p => p.doy.forEach((d,i) => {
      if (d < lo) lo = d; if (d > hi) hi = d;
      [p.ndmi[i], p.ndvi[i]].forEach(v => {
        if (v == null) return;
        if (v < ylo) ylo = v; if (v > yhi) yhi = v; });
    }));
    if (lo < 1e9){ X0 = lo - 4; X1 = hi + 4;
      Y0 = Math.min(0, ylo - 0.03); Y1 = Math.max(1, yhi + 0.03); }
  })();
  let cur = 0;
  let showCuts = true;      // H toggles. Hide it to label blind, then reveal.
  const labels = Object.assign({}, __PRIOR__);   // keep labels from earlier sessions
  DATA.forEach(p => { if (!labels[p.key]) labels[p.key] = p.have.slice(); });

  const svg = document.getElementById('plot');
  const NS = 'http://www.w3.org/2000/svg';
  const sx = d => M.l + (d - X0) / (X1 - X0) * (W - M.l - M.r);
  const xs = v => X0 + (v - M.l) / (W - M.l - M.r) * (X1 - X0);
  const sy = v => H - M.b - (v - Y0) / (Y1 - Y0) * (H - M.t - M.b);

  function el(t, a){ const e = document.createElementNS(NS, t);
    for (const k in a) e.setAttribute(k, a[k]); return e; }

  function draw(){
    const p = DATA[cur];
    svg.innerHTML = '';
    document.getElementById('pos').innerHTML =
      (cur+1) + ' / ' + DATA.length + ' &nbsp; ' + p.pid + ' &nbsp; ' + p.year
      + ' &nbsp; <span style="padding:1px 7px;border-radius:9px;color:#fff;'
      + 'background:' + (p.regime === 'extensive' ? '#8250df' : '#1a7f37')
      + '">LNF ' + p.code + ' &middot; ' + p.regime + '</span>'
      + ' <span style="color:#666;font-weight:400">' + p.expect + '</span>'
      + '<div style="font-family:ui-monospace,monospace;font-size:11px;'
      + 'color:#888;font-weight:400">' + p.key + '</div>';

    const step = (X1 - X0) > 150 ? 20 : 10;
    for (let d = Math.ceil(X0/step)*step; d <= X1; d += step){
      svg.appendChild(el('line', {x1:sx(d), x2:sx(d), y1:M.t, y2:H-M.b,
        stroke:'#eee'}));
      const tx = el('text', {x:sx(d), y:H-M.b+14, 'font-size':10,
        fill:'#666', 'text-anchor':'middle'});
      tx.textContent = d; svg.appendChild(tx);
    }
    [0,0.2,0.4,0.6,0.8,1].filter(v => v >= Y0 && v <= Y1).forEach(v => {
      svg.appendChild(el('line', {x1:M.l, x2:W-M.r, y1:sy(v), y2:sy(v),
        stroke:'#eee'}));
      const t = el('text', {x:M.l-6, y:sy(v)+3, 'font-size':10, fill:'#666',
        'text-anchor':'end'}); t.textContent = v.toFixed(1); svg.appendChild(t);
    });

    // Legal earliest cutting date, extensive meadows only. Anything you label
    // to the left of this line is either an error or a genuine violation.
    if (p.legal){
      svg.appendChild(el('rect', {x:M.l, y:M.t, width:Math.max(0, sx(p.legal)-M.l),
        height:H-M.t-M.b, fill:'#8250df', opacity:0.07}));
      svg.appendChild(el('line', {x1:sx(p.legal), x2:sx(p.legal), y1:M.t,
        y2:H-M.b, stroke:'#8250df', 'stroke-width':1.6}));
      const lt = el('text', {x:sx(p.legal)-5, y:M.t+12, 'font-size':10,
        fill:'#8250df', 'text-anchor':'end'});
      lt.textContent = '15 June'; svg.appendChild(lt);
    }

    // detector proposals
    if (showCuts) p.cuts.forEach(c => svg.appendChild(el('line', {x1:sx(c),
      x2:sx(c), y1:M.t, y2:H-M.b, stroke:'#ff3b30', 'stroke-width':1.4,
      'stroke-dasharray':'5,4', opacity:0.55})));

    // series
    const path = (arr, col) => {
      let dd = '';
      p.doy.forEach((x,i) => { if (!p.haze[i])
        dd += (dd ? 'L' : 'M') + sx(x) + ' ' + sy(arr[i]); });
      svg.appendChild(el('path', {d:dd, fill:'none', stroke:col,
        'stroke-width':1.6}));
      p.doy.forEach((x,i) => {
        if (p.haze[i]) return;
        svg.appendChild(el('circle', {cx:sx(x), cy:sy(arr[i]), r:2.6,
          fill:col}));
      });
    };
    path(p.ndvi, '#9aa5b1');
    path(p.ndmi, '#1f6feb');
    p.doy.forEach((x,i) => { if (p.haze[i])
      svg.appendChild(el('text', {x:sx(x), y:sy(p.ndmi[i])+4, 'font-size':11,
        fill:'#bbb', 'text-anchor':'middle'})).textContent = '\u00d7'; });

    const lg = el('text', {x:M.l+4, y:16, 'font-size':11, fill:'#333'});
    lg.textContent = 'NDMI (blue)  ·  NDVI (grey)  ·  \u00d7 = haze';
    svg.appendChild(lg);

    // user labels
    labels[p.key].forEach(c => {
      svg.appendChild(el('line', {x1:sx(c), x2:sx(c), y1:M.t, y2:H-M.b,
        stroke:'#0b6bcb', 'stroke-width':2}));
      const h = el('rect', {x:sx(c)-5, y:M.t-9, width:10, height:10,
        fill:'#0b6bcb', class:'handle', 'data-d':c, style:'cursor:pointer'});
      svg.appendChild(h);
    });
    dump();
  }

  function nearestObs(d){
    const p = DATA[cur]; let best = p.doy[0], bd = 1e9;
    p.doy.forEach((x,i) => { if (p.haze[i]) return;
      const dd = Math.abs(x - d); if (dd < bd){ bd = dd; best = x; } });
    return best;
  }

  svg.addEventListener('click', e => {
    const p = DATA[cur];
    const r = svg.getBoundingClientRect();
    const px = (e.clientX - r.left) * (W / r.width);
    if (e.target.classList && e.target.classList.contains('handle')){
      const d = +e.target.getAttribute('data-d');
      labels[p.key] = labels[p.key].filter(v => v !== d);
      draw(); return;
    }
    if (px < M.l || px > W - M.r) return;
    const d = nearestObs(Math.round(xs(px)));
    if (!labels[p.key].includes(d)) labels[p.key].push(d);
    labels[p.key].sort((a,b) => a-b);
    draw();
  });

  function dump(){
    // Iterate over ALL keys, not just this session's panels: labels carried in
    // from a previous session must survive into the output, or copying the box
    // would silently discard them.
    const meta = {};
    DATA.forEach(p => meta[p.key] = ' # LNF ' + p.code + ' ' + p.regime);
    let s = 'labels = {\n', n = 0;
    Object.keys(labels).sort().forEach(k => {
      if (!labels[k].length) return;
      s += '    "' + k + '": [' + labels[k].join(', ') + '],'
         + (meta[k] || '') + '\n';
      n += labels[k].length;
    });
    s += '}';
    document.getElementById('out').value = s;
    document.getElementById('cnt').textContent =
      Object.keys(labels).filter(k => labels[k].length).length
      + ' parcel-years, ' + n + ' cuts';
  }

  document.getElementById('prev').onclick = () => {
    cur = (cur - 1 + DATA.length) % DATA.length; draw(); };
  document.getElementById('next').onclick = () => {
    cur = (cur + 1) % DATA.length; draw(); };
  document.getElementById('accept').onclick = () => {
    const p = DATA[cur];
    p.cuts.forEach(c => { const d = nearestObs(Math.round(c));
      if (!labels[p.key].includes(d)) labels[p.key].push(d); });
    labels[p.key].sort((a,b) => a-b); draw(); };
  document.getElementById('clear').onclick = () => {
    labels[DATA[cur].key] = []; draw(); };
  document.getElementById('hide').onclick = () => {
    showCuts = !showCuts;
    document.getElementById('hide').textContent =
      showCuts ? 'hide detector' : 'show detector';
    draw(); };

  document.getElementById('lblwrap').addEventListener('keydown', e => {
    if (e.key === 'ArrowRight') document.getElementById('next').click();
    if (e.key === 'ArrowLeft')  document.getElementById('prev').click();
    if (e.key.toLowerCase() === 'a') document.getElementById('accept').click();
    if (e.key.toLowerCase() === 'c') document.getElementById('clear').click();
    if (e.key.toLowerCase() === 'h') document.getElementById('hide').click();
  });
  document.getElementById('lblwrap').tabIndex = 0;

  draw();
  } catch (e) { fail(e && e.message ? e.message : String(e)); }
})();
</script>
"""
