"""Generate the study-notes badge (assets/logo.svg).

Circular badge (circle = training/learning material) in the teal & cream
palette used for Min-Yao Jhu's learning-notes repositories: a notebook with
coloured chapter tabs and a line of code from the book. Change TITLE and
CODE to reuse it for another book.

Usage: python3 assets/make_logo.py
"""
from pathlib import Path

TITLE = "SeuratExtend"
CODE = ("pbmc |&gt;", "DimPlot2()")   # two lines of code on the front page (XML-escaped)

BG, CREAM, TEAL, TEAL_LT, CORAL, GOLD = "#0F3B3A", "#F3E9D2", "#2E8B83", "#7ED0C4", "#E9967A", "#F2C14E"
FONT = 'font-family="Helvetica, Arial, sans-serif"'
MONO = 'font-family="SFMono-Regular, Menlo, Consolas, monospace"'


def logo(title=TITLE, code=CODE):
    b = [f'<circle cx="120" cy="120" r="114" fill="{BG}" stroke="{CREAM}" stroke-width="7"/>']
    for k, dx in enumerate([10, 5, 0]):                      # stacked notebook pages
        b.append(f'<rect x="{78 + dx}" y="{50 + dx}" width="84" height="96" rx="6" '
                 f'fill="{CREAM if k == 2 else TEAL}" stroke="{BG}" stroke-width="2"/>')
    for i, c in enumerate([CORAL, TEAL_LT, GOLD]):            # chapter tabs
        b.append(f'<rect x="160" y="{62 + i * 18}" width="12" height="12" rx="2" fill="{c}"/>')
    b.append(f'<text x="120" y="96" text-anchor="middle" {MONO} font-size="15" font-weight="bold" fill="{BG}">{code[0]}</text>')
    b.append(f'<text x="120" y="118" text-anchor="middle" {MONO} font-size="13" font-weight="bold" fill="{TEAL}">{code[1]}</text>')
    b.append(f'<text x="120" y="178" text-anchor="middle" {FONT} font-size="23" font-weight="bold" fill="{CREAM}">{title}</text>')
    b.append(f'<text x="120" y="199" text-anchor="middle" {FONT} font-size="12.5" letter-spacing="2.4" fill="{TEAL_LT}">STUDY NOTES</text>')
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 240 240" width="240" height="240" role="img" '
            f'aria-label="{title} study notes badge">\n' + "\n".join(b) + "\n</svg>\n")


if __name__ == "__main__":
    out = Path(__file__).with_name("logo.svg")
    out.write_text(logo())
    print("wrote", out)
