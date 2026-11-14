#!/usr/bin/env python3
"""technote <name>.tex -> <out>.html with pandoc, figures go to <out>/

usage: tex2html.py <name>.tex <out>.html [<pdf link>]
"""
import re, shutil, subprocess, sys
from pathlib import Path

src, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
pdflink = sys.argv[3] if len(sys.argv) > 3 else out.with_suffix(".pdf").name
tex = src.read_text(encoding="utf-8")

# lines starting with '%html ' are active in the html build only
tex = re.sub(r"^%html ", "", tex, flags=re.M)

BR = "ZZBRZZ"   # becomes <br> after pandoc, which drops \newline in cells

# multi-line table cells: \\ inside braces is a line break, not a row end
def cellbreaks(m):
    s, res, depth, i = m[0], [], 0, 0
    while i < len(s):
        if s[i] == "\\":
            if s[i+1] == "\\" and depth > 0:
                i += 2 + (re.match(r"\[.*?\]", s[i+2:]) or re.match("", "")).end()
                res.append(BR)
            else:
                res.append(s[i:i+2])
                i += 2
            continue
        depth += (s[i] == "{") - (s[i] == "}")
        res.append(s[i])
        i += 1
    return "".join(res)
tex = re.sub(r"\\begin\{tabular\}.*?\\end\{tabular\}", cellbreaks, tex, flags=re.S)
tex = re.sub(r"\\makecell(\[.*?\])?|\\raisebox\{[^{}]*\}", "", tex)

# thebibliography -> numbered list, \cite{a,b} -> [1, 2]
keys = re.findall(r"\\bibitem\{(.*?)\}", tex)
tex = re.sub(r"\\begin\{thebibliography\}\{.*?\}", r"\\section*{References}\\begin{enumerate}", tex)
tex = tex.replace(r"\end{thebibliography}", r"\end{enumerate}")
tex = re.sub(r"\\bibitem\{(.*?)\}", lambda m: r"\item \hypertarget{ref-%s}{}" % m[1], tex)

def cite(m):
    refs = []
    for key in (k.strip() for k in m[1].split(",")):
        if key not in keys:
            print(f"{src.name}: unknown citation {key}", file=sys.stderr)
        refs.append(r"\hyperlink{ref-%s}{%s}" % (key, keys.index(key) + 1 if key in keys else "?"))
    return "[" + ", ".join(refs) + "]"
tex = re.sub(r"\\cite\{(.*?)\}", cite, tex)

# caption numbers as in the pdf
count = {}
def caption(m):
    if r"\caption{" not in m[0]:
        return m[0]
    count[m[1]] = count.get(m[1], 0) + 1
    return m[0].replace(r"\caption{", r"\caption{%s %d: " % (m[1].capitalize(), count[m[1]]), 1)
tex = re.sub(r"\\begin\{(table|figure)\*?\}.*?\\end\{\1\*?\}", caption, tex, flags=re.S)
tex = re.sub(r"\\date\{.*", lambda m: m[0].replace(r"\quad", "~"), tex)

# MathJax does not expand text macros inside \text
tex = re.sub(r"\\text\{\\code\{([^{}]*)\}\}", r"\\texttt{\1}", tex)

# figures: browsers need png/svg, take the one next to a pdf figure
figdir = out.with_suffix("")
def figure(m):
    fig = src.parent / m[2]
    if fig.suffix == ".pdf":
        alt = [f for f in (fig.with_suffix(".svg"), fig.with_suffix(".png")) if f.exists()]
        if alt:
            fig = alt[0]
        else:
            print(f"{src.name}: no svg or png next to {fig.name}", file=sys.stderr)
    figdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(fig, figdir / fig.name)
    return r"\includegraphics%s{%s/%s}" % (m[1] or "", figdir.name, fig.name)
tex = re.sub(r"\\includegraphics(\[.*?\])?\{(.*?)\}", figure, tex)

css = """<style>
body { max-width: 52em; }
table { font-size: 0.9em; }
td, th { padding: 0.2em 0.6em; vertical-align: top; }
pre { background: #f6f6f6; padding: 0.6em; border: 1px solid #ddd; }
nav.note { font-size: 0.9em; margin-bottom: 1em; }
span.math.display { display: block; overflow-x: auto; overflow-y: hidden; }
</style>"""
nav = '<nav class="note"><a href="index.html">NGSolve technical notes</a> &middot; <a href="%s">PDF version</a></nav>' % pdflink

# explicit URL: Debian's pandoc defaults to /usr/share/javascript/mathjax
MATHJAX = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml-full.js"
html = subprocess.run(["pandoc", "-f", "latex", "-t", "html5", "-s", "--mathjax=" + MATHJAX, "-V", "lang=en",
                       "-V", "header-includes=" + css, "-V", "include-before=" + nav],
                      input=tex, capture_output=True, text=True, encoding="utf-8", cwd=src.parent)
sys.stderr.write(html.stderr)
if html.returncode:
    sys.exit(html.returncode)
# pandoc 2.x template adds polyfill.io (hijacked domain)
page = re.sub(r'\s*<script src="https://polyfill\.io/[^"]*"></script>', "", html.stdout)
out.write_text(page.replace(BR, "<br>"), encoding="utf-8")
