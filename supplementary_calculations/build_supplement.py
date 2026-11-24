"""Build the accompanying, internally linked mathematical calculations with Pandoc."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parent


def assemble():
    index = (ROOT / "README.md").read_text()
    entries = re.findall(r"^- \[([A-H]\d+)\. .*?\]\(([^)]+\.md)\)$", index, re.M)
    ids = {(ROOT / name).resolve(): label for label, name in entries}
    if len(entries) != len(ids) or set(ids) != {p.resolve() for p in ROOT.glob("*/*.md")}:
        raise ValueError("Each calculation must occur exactly once in the index.")
    _, articles = index.split("## Articles cited\n", 1)
    titles, groups = articles.split("\n## ", 1)
    result = ["# Articles cited", titles.strip()]
    for line in ("## " + groups).splitlines():
        if line.startswith("## "):
            if line == "## Arithmetic":
                break
            result.append(line[1:])
        match = re.match(r"^- \[([A-H]\d+)\. .*?\]\(([^)]+\.md)\)$", line)
        if not match:
            continue
        label, name = match.groups()
        source = ROOT / name
        content = source.read_text().splitlines()
        title = content.pop(0).removeprefix("# ")
        if not title.startswith(label + "."):
            raise ValueError(f"Inconsistent heading in {source}")
        result.append(f"## {title} {{#{label}}}")
        # Preserve the authors' displayed equations, renumbering only heading levels.
        body = re.sub(r"^(#{1,5}) ", r"#\1 ", "\n".join(content), flags=re.M)

        def link(m):
            text, url = m.groups()
            if "://" in url or url.startswith("#"):
                return m.group()
            target = (source.parent / url.partition("#")[0]).resolve()
            if target in ids:
                return f"[{text}](#{ids[target]})"
            if not target.exists():
                raise ValueError(f"Missing local reference {url} in {source}")
            if target.is_relative_to(ROOT):
                return f"[{text}]({target.relative_to(ROOT).as_posix()})"
            return m.group()

        body = re.sub(r"\[([^\]\n]+)\]\(([^)\s]+)\)", link, body)
        body = re.sub(r"\{\\rm\s+([^{}]*)\}", r"\\mathrm{\1}", body)
        body = body.replace(r"\hbox{", r"\text{")
        result.append(body)
    return "\n\n".join(result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pandoc", default="pandoc")
    args = parser.parse_args()
    css = """<style>html {font-size:17px} body {max-width:52rem; margin:3rem auto;
    padding:0 1.5rem;line-height:1.55;color:#18232d} h1,h2{line-height:1.2}
    main,section{min-width:0} math[display=block]{display:block;overflow-x:auto;
    padding:.65rem 0} pre{overflow-x:auto} nav ul{list-style:none}
    a{color:#214b78} @media print{body{max-width:none;margin:0}h1{break-before:page}}
    </style>"""
    with tempfile.TemporaryDirectory() as directory:
        header = Path(directory) / "style.html"
        header.write_text(css)
        subprocess.run(
            [args.pandoc, "--standalone", "--from=markdown", "--to=html5", "--math-method=mathml",
             "--toc", "--toc-depth=2", "--include-in-header", str(header),
             "--metadata", "title=Supplementary calculations for branched covers and spin structures in local curve counting",
             "--metadata", "author=Bernd Johannes Wuebben", "--metadata", "date=24 November 2026",
             "--output", str(ROOT / "supplement.html")],
            input=assemble(), text=True, check=True,
        )
