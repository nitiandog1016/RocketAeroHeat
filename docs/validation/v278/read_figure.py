"""从固定 NASA 原 PDF 重建图 9 两次读数，保留旋转、裁切和倾斜轴标定。"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from PIL import Image, ImageDraw


def calibrate(ticks, symbols, symbol_error):
    # 两个斜轴同时标定。残差与符号位置敏感性只描述读图误差。
    pixels=np.array([row[:2] for row in ticks])
    targets=np.array([row[2:] for row in ticks])
    matrix=np.column_stack((pixels,np.ones(len(pixels))))
    transform=np.linalg.lstsq(matrix,targets,rcond=None)[0]
    residual=abs(matrix@transform-targets).max(axis=0)
    points=np.column_stack((symbols,np.ones(len(symbols))))@transform
    spread=residual[1]+np.array(symbol_error)*np.abs(transform[:2,1]).sum()
    return points,spread,dict(ticks=ticks,symbols=symbols,symbolHalfWidthPixels=symbol_error,
                              affinePixelToAngleNu=transform.tolist(),maxTickResidual=residual.tolist())


def digitize(pdf,out):
    digest=hashlib.sha256(pdf.read_bytes()).hexdigest()
    assert digest=="48f470f84a8d7ed11e9c3f4739f98ef119db5c39eec41869e6fde693e13967ea"
    # 低倍读图直接取菱形中心。高倍用上下顶点中点，零度与纵轴和理论线重叠。
    primary=[(121,417),(201,413),(281,410),(361,409),(442,406),(523,409)]
    tips=[(122,626,666),(241,619,659),(359,614,654),(482,613,654),(601,608,649),(726,612,653)]
    secondary=[(x,.5*(top+bottom)) for x,top,bottom in tips]
    ticks1=[(126,537,0,0),(446.5,522,40,0),(769.5,508.3,80,0),
            (116,283.6,0,800),(110.6,157.4,0,1200),(105.4,29,0,1600)]
    ticks2=[(368.5,814,20,0),(609.5,803,40,0),(850,792,60,0),
            (122.5,637,0,400),(114.5,444,0,800),(106.5,255.5,0,1200)]
    out.mkdir(parents=True,exist_ok=False)
    readings=[]
    metadata=[]
    for dpi,box,ticks,points,err in [
        (240,(1500,940,2320,1550),ticks1,primary,[7,4,4,4,5,4]),
        (360,(2310,1390,3420,2240),ticks2,secondary,[10,5,5,5,7,5])]:
        prefix=out/str(dpi)
        subprocess.run(["pdftoppm","-f","14","-l","14","-singlefile","-gray","-r",str(dpi),"-png",str(pdf),str(prefix)],check=True)
        upright=Image.open(prefix.with_suffix(".png")).rotate(270,expand=True)
        if dpi==240:
            upright.crop((1430,690,2360,1680)).save(out/"source_figure9.png")
        image=upright.crop(box).convert("RGB")
        draw=ImageDraw.Draw(image)
        for x,y,*_ in ticks:
            draw.ellipse((x-5,y-5,x+5,y+5),outline="#236383",width=2)
        for x,y in points:
            draw.ellipse((x-15,y-15,x+15,y+15),outline="#b13f1f",width=2)
        image.save(out/("markers_%s.png"%dpi))
        values,spread,info=calibrate(ticks,points,err)
        info.update(dpi=dpi,cropAfterClockwiseRotation=box)
        readings.append((values,spread))
        metadata.append(info)
    (a,ea),(b,eb)=readings
    nominal=np.arange(0,51,10)
    assert np.max(abs(a[:,0]-nominal))<1 and np.max(abs(b[:,0]-nominal))<1
    assert np.max(abs(a[:,1]-b[:,1]))<10
    # 表 III 提供的是另一种数据整理方式。它与测得菱形及独立理论计算分列。
    fit_a=[1.2436,.9734,.8373,.8160,.6569,.6460]
    fit_b=[.4774,.5008,.5125,.5109,.5241,.5218]
    with (out/"stations.csv").open("x",encoding="ascii",newline="") as stream:
        writer=csv.writer(stream,lineterminator="\n")
        writer.writerow(["angleDeg","figureNu","readLowNu","readHighNu","secondReadNu","tableA","tableB"])
        for i,angle in enumerate(nominal):
            writer.writerow([angle,*[round(v,4) for v in (a[i,1],min(a[i,1]-ea[i],b[i,1]-eb[i]),
                             max(a[i,1]+ea[i],b[i,1]+eb[i]),b[i,1])],fit_a[i],fit_b[i]])
    (out/"calibration.json").write_text(json.dumps(dict(pdfSHA256=digest,figure=9,pdfPage=14,
        highResolutionDiamondTips=tips,readings=metadata),indent=2),encoding="utf8")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf",type=Path)
    parser.add_argument("new_output",type=Path)
    args=parser.parse_args()
    digitize(args.pdf,args.new_output)
