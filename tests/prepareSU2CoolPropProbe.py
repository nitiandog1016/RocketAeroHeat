"""准备一次固定官方运行包的CoolProp能力探测；不使用dryrun或改变箭体配对物理。"""
from pathlib import Path
import hashlib,json,os,re,uuid
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=ROOT/'outputs/v105_validation/plate03';folder=ROOT/'outputs/v120_validation/coolprop_probe';assert not folder.exists()
    # 固定源码要求可压缩非理想气体用ROE/HLLC/JST，避免AUSM配置拒绝掩盖真正能力。
    configSource=ROOT/'tools/runtime/su2-thermal-build/source/Common/src/CConfig.cpp'
    fluidSource=ROOT/'tools/runtime/su2-thermal-build/source/SU2_CFD/src/fluid/CCoolProp.cpp'
    text=(source/'case.cfg').read_text()
    for key,value in [('FLUID_MODEL','COOLPROP'),('CONV_NUM_METHOD_FLOW','ROE'),('ITER','1'),('OUTPUT_WRT_FREQ','1')]:
        text,count=re.subn(r'^'+key+r'=.*$',key+'= '+value,text,flags=re.M);assert count==1
    for key in ['MARKER_FAR','MARKER_SUPERSONIC_OUTLET']:
        text,count=re.subn(r'^'+key+r'=.*\n','',text,flags=re.M);assert count==1
    text+='\nFLUID_NAME= Air\nMARKER_RIEMANN= (farfield, STATIC_SUPERSONIC_INFLOW_PT, 5000, 300, -4, 0, 0, outlet, STATIC_PRESSURE, 5000, 0, 0, 0, 0)\n'
    folder.mkdir();(folder/'case.cfg').write_text(text,encoding='utf-8')
    os.link(str(source/'mesh.su2'),str(folder/'mesh.su2'))
    manifest=dict(kind='fixed-runtime-compressible-property-capability-probe',requiredRuntimeFlavor='official',
                  purpose='test CoolProp instantiation and at most one small real iteration; no physical qualification',
                  sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),configSource,fluidSource,source/'case.cfg',source/'mesh.su2']})
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=60,dryRun=False,
                 jobToken=str(uuid.uuid4()),sourceFiles=[str(p.relative_to(ROOT)) for p in [Path(__file__).resolve(),folder/'case_manifest.json',configSource,fluidSource]])
    (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n',encoding='utf-8');print(str(folder))
if __name__=='__main__':main()
