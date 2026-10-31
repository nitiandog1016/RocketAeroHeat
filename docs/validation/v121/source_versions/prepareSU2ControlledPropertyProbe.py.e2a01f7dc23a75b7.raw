"""复用v120已核实的物性能力算例，改用明确声明的受控运行包；无物理资格。"""
from pathlib import Path
import json,hashlib,os,uuid
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    source=ROOT/'outputs/v120_validation/coolprop_probe';folder=ROOT/'outputs/v121_validation/capability01'
    assert not folder.exists();folder.mkdir()
    for name in ['case.cfg','mesh.su2']:os.link(str(source/name),str(folder/name))
    manifest=dict(kind='controlled-compressible-property-capability-probe',requiredRuntimeFlavor='rah-coolprop-v1',
       thermodynamicModel='CoolProp-7.2.0-HEOS-Air',physicsQualified=False,
       sourceSHA256={str(p):sha(p) for p in [Path(__file__).resolve(),source/'case.cfg',source/'mesh.su2',source/'engine_process.json']})
    (folder/'case_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    request=dict(workDirectory=str(folder),configFile='case.cfg',meshFile='mesh.su2',threads=2,maximumSeconds=60,dryRun=False,jobToken=str(uuid.uuid4()),
      sourceFiles=[str(p.relative_to(ROOT)) for p in [Path(__file__).resolve(),folder/'case_manifest.json',ROOT/'outputs/v121_validation/runtime_installed.json']])
    (folder/'engine_request.json').write_text(json.dumps(request,indent=2)+'\n');print(str(folder))
if __name__=='__main__':main()
