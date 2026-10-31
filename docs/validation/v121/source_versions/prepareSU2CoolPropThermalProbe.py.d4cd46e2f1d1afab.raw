"""第0步零热流不能验证导热；准备有预算的同网格20步非零热流校核。"""
from pathlib import Path
import hashlib,json,os,re,uuid
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=ROOT/'outputs/v121_validation/capability01';dest=ROOT/'outputs/v121_validation/thermal20'
    assert not dest.exists();dest.mkdir()
    cfg=(source/'case.cfg').read_text()
    for key,value in [('ITER','20'),('OUTPUT_WRT_FREQ','10'),('WRT_RESTART_OVERWRITE','NO'),('WRT_SURFACE_OVERWRITE','NO')]:
        cfg,n=re.subn('^'+key+r'=.*$',key+'= '+value,cfg,flags=re.M);assert n==1
    (dest/'case.cfg').write_text(cfg);os.link(str(source/'mesh.su2'),str(dest/'mesh.su2'))
    manifest=dict(kind='nonzero-native-coolprop-wall-gradient-check',requiredRuntimeFlavor='rah-coolprop-v1',
        thermodynamicModel='CoolProp-7.2.0-HEOS-Air',physicsQualified=False,
        originalZeroFluxProbe=str(source),purpose='20 startup steps test nonzero wall-gradient data path; no physical qualification',
        sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),source/'case.cfg',source/'mesh.su2',source/'wall_audit.json',ROOT/'docs/validation/v121/PLAN.md']})
    (dest/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    request=dict(workDirectory=str(dest),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=60,dryRun=False,jobToken=str(uuid.uuid4()),
        sourceFiles=[str(p.relative_to(ROOT)) for p in [Path(__file__).resolve(),dest/'case_manifest.json',ROOT/'docs/validation/v121/PLAN.md']])
    (dest/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n');print(str(dest))
if __name__=='__main__':main()
