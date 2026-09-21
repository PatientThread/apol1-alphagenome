#!/usr/bin/env python3
"""
format_for_submission.py

Applies the journal's general submission formatting to MANUSCRIPT.docx, which
the markdown converter does not do: double spacing, page numbers, continuous
line numbers, a plain unshaded table without thousands separators, and the
table legend placed below the table.

Run AFTER md_to_docx.py. Writes MANUSCRIPT_submission.docx and leaves the
converter output untouched.

Author: Christopher Lawrence
"""

from __future__ import annotations

import copy
import re
from pathlib import Path

from docx import Document
from docx.enum.text import WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt

BODY_FONT = "Times New Roman"   # universally available; Calibri is a Windows default
CODE_FONT = "Courier New"       # Menlo is macOS-only and substitutes unpredictably
BODY_SIZE = Pt(12)

HERE = Path(__file__).resolve().parent
SRC = HERE / "MANUSCRIPT.docx"
OUT = HERE / "MANUSCRIPT_submission.docx"


def add_line_numbers(section) -> None:
    """Continuous line numbering, restarting at each section."""
    sectPr = section._sectPr
    ln = sectPr.find(qn("w:lnNumType"))
    if ln is None:
        ln = OxmlElement("w:lnNumType")
        sectPr.append(ln)
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:start"), "1")
    ln.set(qn("w:restart"), "continuous")


def add_page_numbers(section) -> None:
    """PAGE field centred in the footer."""
    p = section.footer.paragraphs[0] if section.footer.paragraphs \
        else section.footer.add_paragraph()
    p.alignment = 1
    for text, typ in (("begin", "w:fldChar"), ("PAGE", "w:instrText"),
                      ("end", "w:fldChar")):
        r = OxmlElement("w:r")
        el = OxmlElement(typ)
        if typ == "w:fldChar":
            el.set(qn("w:fldCharType"), text)
        else:
            el.set(qn("xml:space"), "preserve")
            el.text = " PAGE "
        r.append(el)
        p._p.append(r)


def strip_shading(table) -> None:
    """Remove fills and borders styling; leave a plain grid."""
    table.style = "Table Grid"
    for row in table.rows:
        for cell in row.cells:
            tcPr = cell._tc.get_or_add_tcPr()
            for shd in tcPr.findall(qn("w:shd")):
                tcPr.remove(shd)
            shd = OxmlElement("w:shd")
            shd.set(qn("w:val"), "clear")
            shd.set(qn("w:fill"), "auto")
            tcPr.append(shd)
            for para in cell.paragraphs:
                for run in para.runs:
                    run.font.color.rgb = None


def repeat_header(table) -> None:
    trPr = table.rows[0]._tr.get_or_add_trPr()
    if trPr.find(qn("w:tblHeader")) is None:
        h = OxmlElement("w:tblHeader")
        h.set(qn("w:val"), "true")
        trPr.append(h)


def strip_thousands(table) -> None:
    """29,989 -> 29989 in table cells only; prose keeps separators."""
    pat = re.compile(r"(?<=\d),(?=\d{3}\b)")
    for row in table.rows:
        for cell in row.cells:
            for para in cell.paragraphs:
                for run in para.runs:
                    if "," in run.text:
                        run.text = pat.sub("", run.text)


def main() -> None:
    doc = Document(SRC)

    # One body font throughout, and one monospace font for code identifiers.
    # The converter emits Menlo for markdown backticks, which exists only on
    # macOS; an editor opening this on Windows would get an arbitrary
    # substitution mid-sentence.
    normal = doc.styles["Normal"]
    normal.font.name = BODY_FONT
    normal.font.size = BODY_SIZE
    rpr = normal.element.get_or_add_rPr()
    rf = rpr.find(qn("w:rFonts"))
    if rf is None:
        rf = OxmlElement("w:rFonts")
        rpr.append(rf)
    for a in ("w:ascii", "w:hAnsi", "w:cs"):
        rf.set(qn(a), BODY_FONT)

    def restyle(runs):
        for run in runs:
            name = run.font.name
            if name is None:
                r = run._r.find(qn("w:rPr"))
                if r is not None:
                    f = r.find(qn("w:rFonts"))
                    if f is not None:
                        name = f.get(qn("w:ascii")) or f.get(qn("w:hAnsi"))
            mono = name in ("Menlo", "Consolas", "Courier", "Monaco",
                            "DejaVu Sans Mono")
            run.font.name = CODE_FONT if mono else BODY_FONT
            run.font.size = BODY_SIZE
            r = run._r.get_or_add_rPr()
            f = r.find(qn("w:rFonts"))
            if f is None:
                f = OxmlElement("w:rFonts")
                r.append(f)
            for a in ("w:ascii", "w:hAnsi", "w:cs"):
                f.set(qn(a), CODE_FONT if mono else BODY_FONT)

    for para in doc.paragraphs:
        pf = para.paragraph_format
        pf.line_spacing_rule = WD_LINE_SPACING.DOUBLE
        pf.space_after = Pt(0)
        restyle(para.runs)

    for section in doc.sections:
        add_line_numbers(section)
        add_page_numbers(section)

    for table in doc.tables:
        strip_shading(table)
        repeat_header(table)
        strip_thousands(table)
        for row in table.rows:
            for cell in row.cells:
                for para in cell.paragraphs:
                    para.paragraph_format.line_spacing_rule = \
                        WD_LINE_SPACING.SINGLE
                    restyle(para.runs)

    # Move the Table 1 caption below the table: the converter emits it above.
    body = doc.element.body
    caps = [p for p in doc.paragraphs
            if p.text.strip().startswith("Table 1.")]
    if caps and doc.tables:
        cap = caps[0]._p
        tbl = doc.tables[0]._tbl
        body.remove(cap)
        tbl.addnext(cap)

    doc.save(OUT)
    print(f"wrote {OUT.name}")
    print(f"  paragraphs double spaced: {len(doc.paragraphs)}")
    print(f"  tables plainified: {len(doc.tables)}")
    print("  line numbers: continuous, page numbers: footer")
    print(f"  fonts: {BODY_FONT} 12pt body, {CODE_FONT} for code identifiers")


if __name__ == "__main__":
    main()
