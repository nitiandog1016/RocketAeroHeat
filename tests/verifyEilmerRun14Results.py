"""不调用主回读器：用Decimal逐行核对原生载荷、CSV、测站内插和全周热功率。"""
from pathlib import Path
from decimal import Decimal,localcontext
import argparse
import csv
import hashlib
import json
import math
import time


def verify(case,result,out):
    started=time.monotonic();case=Path(case);result=Path(result);out=Path(out)
    assert not out.exists()
    meta=[line.split() for line in (case/'lmrsim/loads/loads-metadata').read_text().splitlines()[1:] if line.strip()]
    native=[]
    for block in (0,1):
        files=list((case/'lmrsim/loads'/meta[-1][0]).glob('blk-%04d-*-wall.dat'%block))
        assert len(files)==1
        lines=[line.split() for line in files[0].read_text().splitlines() if line.strip()]
        native.extend(dict(zip(lines[0],row)) for row in lines[1:])
    wall=list(csv.DictReader((result/'wall.csv').open()))
    stations=list(csv.DictReader((result/'stations.csv').open()))
    report=json.loads((result/'report.json').read_text())
    assert len(native)==len(wall)==320 and report['physicsQualified'] is False
    with localcontext() as ctx:
        ctx.prec=45
        D=Decimal;total=D(0);area=D(0);max_error=D(0)
        for n,w in zip(native,wall):
            assert abs(D(n['pos.x'])-D(w['x_m']))<D('3e-14')
            expected=-D(n['outsign'])*(D(n['q_cond'])+D(n['q_diff']))
            error=abs(expected-D(w['heat_into_wall_W_m2']))/max(D(1),abs(expected))
            assert error<D('2e-13');max_error=max(max_error,error)
            assert abs(D(n['p'])-D(w['pressure_Pa']))<D('2e-10')
            assert D(w['Tw_given_K'])==D('295.2')
            # 原生轴对称面积为每弧度；完整周向必须乘2π。
            assert abs(D(w['full_ring_area_m2'])/(D(n['area'])*D(str(2*math.pi)))-1)<D('3e-15')
            total+=D(w['full_ring_area_m2'])*expected;area+=D(w['full_ring_area_m2'])
        assert abs(total-D(str(report['researchHeatPowerIntoWallW'])))/abs(total)<D('3e-15')
        assert abs(area-D(str(report['fullRingWallAreaM2'])))<D('3e-16')
        matched=0;unknown=[];maximum_interpolation_error=D(0)
        for s in stations:
            if not s['prediction']:
                unknown.append(s['sensor']);continue
            x=D(s['x_m']);block=0 if x<D('.1017') else 1
            rows=[w for w in wall if int(float(w['block']))==block]
            found=False
            for a,b in zip(rows[:-1],rows[1:]):
                xa,xb=D(a['x_m']),D(b['x_m'])
                if xa<=x<=xb:
                    key='heat_into_wall_W_m2' if s['quantity']=='q' else 'pressure_Pa'
                    y=D(a[key])+(D(b[key])-D(a[key]))*(x-xa)/(xb-xa)
                    err=abs(y-D(s['prediction']))/max(D(1),abs(y))
                    assert err<D('5e-14');maximum_interpolation_error=max(maximum_interpolation_error,err)
                    found=True;break
            assert found
            matched+=1
        assert matched==76 and unknown==['P30','P32']
    result=dict(passed=True,nativeFaces=320,decimalStationMatches=matched,unmapped=unknown,
        maxNativeHeatRelativeDifference=float(max_error),maxDecimalInterpolationRelativeDifference=float(maximum_interpolation_error),
        decimalHeatPowerW=str(total),decimalAreaM2=str(area),physicsQualified=False,newQualifiedAreaM2=0,
        comparedReportSHA256=hashlib.sha256((result/'report.json').read_bytes()).hexdigest(),
        elapsedSeconds=time.monotonic()-started)
    out.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--result',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();print(json.dumps(verify(a.case,a.result,a.out)))
