"""Figure 3 as a vector composite. Panels A-C (Figure_3_ABC_restyled.pdf, from
Figure_3_and_S3_restyled.R) above panel D (Figure_3D_logos.pdf, from Figure_3D_logos.R on the final importances); both are cairo
vector PDFs of the same width (878 pt) and are stacked page on page with pypdf (no rasterisation). The PNG is the composite PDF
rendered at 300 dpi with pdftoppm.
Writes outputs/figures/Figure_3.{pdf,png} (backups first).
usage: python3 4_visualisation/figures/Fig3_S3/assemble_figure3.py
"""
import os
import subprocess
import sys

from pypdf import PdfReader, PdfWriter, Transformation

HERE = os.path.dirname(os.path.abspath(__file__))
R2 = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(R2, "3_analysis", "integration"))
import w2  # noqa: E402

OUT = os.path.join(R2, "outputs/figures")
ABC, D = os.path.join(HERE, "Figure_3_ABC_restyled.pdf"), os.path.join(HERE, "Figure_3D_logos.pdf")
GAP = 10   # pt between C and D
pa, pd_ = PdfReader(ABC).pages[0], PdfReader(D).pages[0]
wa, ha = float(pa.mediabox.width), float(pa.mediabox.height)
wd, hd = float(pd_.mediabox.width), float(pd_.mediabox.height)
if abs(wa - wd) > 1:
    raise SystemExit("CHECK FAILED page widths differ: A-C %.1f pt, D %.1f pt" % (wa, wd))
w = PdfWriter()
page = w.add_blank_page(width=wa, height=ha + GAP + hd)
page.merge_transformed_page(pd_, Transformation().translate(0, 0))
page.merge_transformed_page(pa, Transformation().translate(0, hd + GAP))
pdf, png = os.path.join(OUT, "Figure_3.pdf"), os.path.join(OUT, "Figure_3.png")
for f in (pdf, png):
    w2.backup(f)
with open(pdf, "wb") as fh:
    w.write(fh)
subprocess.run([os.environ.get("PDFTOPPM", "pdftoppm"), "-png", "-r", "300", "-singlefile", pdf, png[:-4]], check=True)
imgs = subprocess.run([os.environ.get("PDFIMAGES", "pdfimages"), "-list", pdf], capture_output=True, text=True).stdout.strip().splitlines()[2:]
print("Figure 3: %.0f x %.0f pt, raster images in the PDF: %d" % (wa, ha + GAP + hd, len(imgs)))
w2.provenance(OUT, "4_visualisation/figures/Fig3_S3/assemble_figure3.py", [pdf, png], [ABC, D])
