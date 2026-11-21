"""重取论文第 815 页，并按原图符号与刻度复算包围 6 K 的温差区间。"""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess

from PIL import Image, ImageChops

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT / "outputs/v282_validation/source01"
PDF_SHA = "41e9098cdfafea78a5a613a1e81c165d1cb910c6075f106b9fede8fb78b94976"


def read_point(xy, ticks, logarithmic):
    """分别扰动点中心和两个刻度，保留手读符号及线宽引起的横坐标范围。"""
    def temperature(x, left, right):
        if logarithmic:
            return 10**((x-left)/(right-left))
        return 200*(x-left)/(right-left)

    x = xy[0]
    left, right = ticks
    limits = [temperature(x+d, left+a, right+b) for d,a,b in
              product((-2,2), (-.75,.75), (-.75,.75))]
    return dict(pixelXY=xy, deltaTK=temperature(x,left,right), readRangeK=[min(limits),max(limits)])


def run(pdf, out):
    assert hashlib.sha256(pdf.read_bytes()).hexdigest()==PDF_SHA
    out.mkdir(parents=True,exist_ok=False)
    subprocess.run(["pdftoppm","-f","6","-l","6","-r","230","-png","-singlefile",
                    str(pdf),str(out/"page815")],check=True,timeout=30)
    image=Image.open(out/"page815.png").convert("RGB")
    assert image.size==(1901,2690)
    prior=json.loads((HERE.parent/"v282/evidence/comparison.json").read_text(encoding="utf8"))["experiment"]
    graph=SOURCE/"figure7.png"
    assert hashlib.sha256(graph.read_bytes()).hexdigest()==prior["cropSHA256"]
    figures=[]
    for number,box,ticks,xy in (
        (7,[1064,1479,1691,1936],prior["xPixelsAt0and200K"],[(101.5,263.5),prior["pointXY"]]),
        (8,[975,1985,1740,2410],[264.5,413.5],[(337,152),(425,151)])):
        crop=image.crop(box)
        crop.save(out/f"figure{number}.png")
        if number==7:
            assert ImageChops.difference(crop,Image.open(graph).convert("RGB")).getbbox() is None
        points=[read_point(p,ticks,number==8) for p in xy]
        figures.append(dict(figure=number,boxOnPage=box,xTicksPx=ticks,xTicksK=[1,10] if number==8 else [0,200],
                            axis="log10" if number==8 else "linear",points=points))
    assert figures[0]["points"][1]["readRangeK"]==prior["deltaTReadRangeK"]
    gap=[max(f["points"][0]["readRangeK"][1] for f in figures),
         min(f["points"][1]["readRangeK"][0] for f in figures)]
    assert gap[0]<6<gap[1]
    record=dict(doi=prior["doi"],pdfSHA256=PDF_SHA,printedPage=815,specimen="B",speedMps=3,
                nominalBeadAreaM2=3.094e-6,pointHalfWidthPx=2,tickHalfWidthPx=.75,
                dpi=230,figures=figures,openGapK=gap,reportedPointUncertaintyWm2K=None)
    (out/"beadExperimentGap.json").write_bytes((json.dumps(record,indent=2)+"\n").encode("utf8"))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    parser.add_argument("--pdf",type=Path,default=SOURCE/"tsubouchi1959.pdf")
    args=parser.parse_args()
    run(args.pdf,args.output)
