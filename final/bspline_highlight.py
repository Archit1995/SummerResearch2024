import re, html
from pathlib import Path

RULES = [
    (r'In [A-Z]+', 'section'),
    (r'sigma_\[\d+\]', 'sigma'),
    (r'Roots.*', 'roots'),
    (r'Success', 'success'),
    (r'failed', 'fail'),
    (r'quotient|remainder|counter', 'step'),
    (r'={5,}|_{5,}', 'divider'),
    (r'[0-9]','purpletext'),
    (r'checking', 'cyantext'),   # <-- NEW: example pattern
]

def highlight_line(line):
    esc = html.escape(line)
    for pat, cls in RULES:
        esc = re.sub(pat, lambda m: f'<span class="{cls}">{m.group(0)}</span>', esc)
    return esc

def main(inp, outp):
    lines = Path(inp).read_text(errors="ignore").splitlines()
    out = ['<html><head><style>',
           'body{background:#0d1117;color:#e6edf3;font-family:monospace;white-space:pre;}',
           '.section{color:#58a6ff;font-weight:bold}',
           '.sigma{color:#ffcc00}',
           '.roots{color:#c3e88d}',
           '.success{color:#00ff00;font-weight:bold}',
           '.fail{color:#ff5555;font-weight:bold}',
           '.step{color:#ff79c6}',
           '.divider{color:#888}',
           '.lineno{color:#555;margin-right:1em;}',
           '.cyantext{color:#00ffff;font-weight:bold}',   # <-- NEW STYLE
           '.purpletext{color:green;font-weight:bold}',
           '</style></head><body>']

    for i, l in enumerate(lines, start=1):
        numbered = f'<span class="lineno">{i:6d}:</span> {highlight_line(l)}'
        out.append(numbered)

    out.append('</body></html>')
    Path(outp).write_text('\n'.join(out), encoding='utf-8')

# Run on your file
main("bspline3.txt", "bspline3_highlight.html")

