"""从原始运行14列独立重算测站/单位；不调用准备程序及其换算函数。"""
from pathlib import Path
import argparse,csv,json,hashlib,math
from decimal import Decimal,localcontext

def verify(folder,output):
    root=Path(__file__).resolve().parents[1];raw=root/'tests/data/eilmer_cylinder_flare/raw'
    assert not output.exists();output.mkdir(parents=True)
    files=[folder/'heat_flux_SI.csv',folder/'pressure_SI.csv',folder/'conditions.json']
    files += [raw/'cylinder-extended-flare-heat-transfer.data',raw/'cylinder-extended-flare-pressure.data']
    before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    report={'passed':False,'nativeTasks':0,'newQualifiedArea_m2':0}
    try:
        with localcontext() as ctx:
            ctx.prec=50
            # 独立单位定义：国际表Btu、英尺、磅力和英寸均直接以SI定义。
            factors=[Decimal('1055.05585262')/Decimal('0.09290304'),
                     Decimal('4.4482216152605')/Decimal('0.00064516')]
            rows_total=0;max_error=0;by_sensor={}
            for raw_path,csv_path,factor in zip(files[3:],files[:2],factors):
                expected=[]
                for line in raw_path.read_text(encoding='utf8').splitlines():
                    tokens=line.lstrip('# ').split()
                    if tokens and ((tokens[0].startswith('HT') and tokens[0][2:].isdigit()) or
                                   (tokens[0].startswith('P') and tokens[0][1:].isdigit())):
                        assert len(tokens)==10
                        expected.append((tokens,line.lstrip().startswith('#')))
                with csv_path.open(encoding='utf-8-sig',newline='') as stream:actual=list(csv.DictReader(stream))
                assert len(actual)==len(expected)
                for (tokens,comment),row in zip(expected,actual):
                    rows_total+=1;by_sensor[row['sensor']]=row
                    assert row['sensor']==tokens[0] and row['source_value']==tokens[-1]
                    x=Decimal(tokens[1])*Decimal('0.1017016')
                    assert Decimal(row['x_from_leading_edge_m'])==x and Decimal(row['project_X_m'])==-x
                    assert (row['record_status']=='source-commented')==comment
                    assert row['physical_qualification']=='not-assessed' and row['uncertainty_SI']==''
                    assert (row['geometry_region']=='outside')==(x>Decimal('.220'))
                    if 'x' in tokens[-1].lower():assert row['value_SI']==''
                    else:
                        truth=Decimal(tokens[-1])*factor;error=float(abs(Decimal(row['value_SI'])-truth)/max(Decimal(1),abs(truth)))
                        assert error<1e-35;max_error=max(max_error,error)
            assert rows_total==81
            assert by_sensor['HT36']['source_value']=='21.01'
            assert by_sensor['P1']['source_value']=='0.008319'
            assert {name for name,row in by_sensor.items() if row['geometry_region']=='outside' and row['record_status']=='source-active'}=={'P30','P32'}
        conditions=json.loads(files[2].read_text(encoding='utf8'))
        d=conditions['codeDerived'];r=conditions['reported'];mu=d['viscosity_Pas']
        assert r['density_kg_m3']==8.810e-4 and r['vibrationTemperature_K']==2467
        assert abs(d['density_kg_m3']*d['gasConstant_J_kgK']*r['temperature_K']/r['pressure_Pa']-1)<1e-14
        assert math.isclose(d['Re_L_originalTable'],d['density_kg_m3']*2304*.1017016/mu,rel_tol=1e-14)
        assert conditions['transport']['actual']['actualEvaluationTemperature_K']==200
        assert not conditions['experimentalQualification'] and conditions['newQualifiedArea_m2']==0
        report.update(passed=True,allOriginalRows=rows_total,maxRelativeUnitError=max_error,
                      peakHeatFlux_W_m2=float(by_sensor['HT36']['value_SI']),
                      outsideActivePreserved=['P30','P32'],commentsAndMissingPreserved=True)
    finally:
        report['inputsStable']=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in before.items())
        (output/'report.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    assert report['passed'] and report['inputsStable'];print(json.dumps(report))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();verify(a.input,a.output)
