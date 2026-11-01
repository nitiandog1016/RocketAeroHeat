"""归档v135准确实现版本、小型导入和真实失败；不删除数据或重跑数值。"""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess,zipfile
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'outputs/v135_validation';DEST=ROOT/'docs/validation/v135'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,d):
    assert not p.exists(),str(p);p.write_bytes((json.dumps(d,ensure_ascii=False,indent=2)+'\n').encode('utf8'))

def old_versions():
    # 运行后的纯说明修改逆向恢复；必须与执行回执SHA完全一致，不靠猜测追认。
    receipt=read(BASE/'v135_interface01_process.json');destination=DEST/'source_versions';destination.mkdir(exist_ok=True)
    replacements={
        'readSU2HybridVolume.m':('棱柱面心积分、金字塔全参考Jacobian和四面体行列式核对完整域。','棱柱中心扇形面积积分与四面体行列式独立核对体积和完整域。'),
        'scripts/su2SSTInterface.py':('棱柱/金字塔/四面体/六面体和二维四边形支持。','棱柱/四面体/六面体和二维四边形支持。')}
    for name in list(replacements)+['rocketAeroHeatWorkbench.m']:
        expected=receipt['sourceSHA256'][name];out=destination/(Path(name).name+'.'+expected+'.raw')
        if out.exists():assert sha(out)==expected;continue
        b=(ROOT/name).read_bytes() if name in replacements else subprocess.check_output(['git','show','0710a6837d42d0595db49a8d43ce5f473140f886:'+name],cwd=ROOT)
        if name in replacements:new,old=replacements[name];b=b.replace(new.encode(),old.encode())
        lf=b.replace(b'\r\n',b'\n');variants=[b,lf,lf.replace(b'\n',b'\r\n')]
        correct=next(v for v in variants if hashlib.sha256(v).hexdigest()==expected);out.write_bytes(correct)
        print('Preserved execution version: '+name)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');args=parser.parse_args();old_versions()
    if args.prepare:return
    assert not (DEST/'delivery_manifest.json').exists()
    inventory=[];unique={}
    for p in sorted(BASE.rglob('*')):
        if not p.is_file():continue
        st=p.stat();assert st.st_ino;unique[(st.st_dev,st.st_ino)]=st.st_size
        inventory.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=st.st_size,sha256=sha(p),fileID=st.st_ino,device=st.st_dev))
    write(DEST/'storage_inventory.json',dict(scope='v135 only',files=inventory,pathCount=len(inventory),distinctFileIDs=len(unique),
        pathLengthBytes=sum(v['bytes'] for v in inventory),distinctFileLengthBytes=sum(unique.values()),isAllocatedDiskUsage=False,
        deletedFiles=0,historicalOutputsTouched=False))
    copies=[]
    def copy(p,target):
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():assert sha(p)==sha(target)
        else:shutil.copyfile(p,target)
        assert sha(p)==sha(target);copies.append(dict(source=p.relative_to(ROOT).as_posix(),archived=target.relative_to(ROOT).as_posix(),bytes=target.stat().st_size,sha256=sha(target)))
    for row in inventory:
        p=ROOT/row['path']
        if row['bytes']<2*1024**2 and p.suffix in ['.json','.log','.mat','.cfg','.su2','.csv','.vtk']:copy(p,DEST/'evidence'/p.relative_to(BASE))
    names=['scripts/meshSU2HybridTransition.py','scripts/su2PyramidGeometry.py','scripts/auditSU2HybridTopology.py',
        'scripts/analyzeSU2DualVolumes.py','scripts/su2SSTInterface.py','scripts/reviewSU2PyramidTransition.py',
        'scripts/meshSU2HybridVolume.py','scripts/optimizeSU2FanExterior.py','scripts/prepareSU2CheckpointContinuation.py',
        'scripts/runMatlabResearch.py','scripts/archiveSU2V135Evidence.py','scripts/runEmbeddedSU2.py',
        'readSU2HybridVolume.m','writeSU2RocketCase.m','su2PyramidGeometry.m','runEmbeddedSU2Case.m',
        'buildSU2WallProjection.m','mapSU2WallResult.m','rocketAeroHeatWorkbench.m','buildSU2CreaseFanSlice.m',
        'buildSU2MeshQualityPreview.m','tests/testSU2PyramidGeometry.py','tests/testSU2PyramidInterface.m',
        'tests/verifySU2PyramidInterface.py','tests/testSU2PyramidWorkbench.m']
    sources=[]
    for name in names:
        p=ROOT/name;target=DEST/'source_versions'/(p.name+'.'+sha(p)+'.raw');copy(p,target)
        sources.append(dict(path=name,sha256=sha(p),archived=target.relative_to(ROOT).as_posix()))
    tasks=[]
    for path in sorted(BASE.glob('*_process.json')):
        d=read(path)
        if 'sourceSHA256' not in d:continue
        assert d['status']=='terminal' and all(d['sourceStableAfter'].values())
        differences=[]
        for name,h in d['sourceSHA256'].items():
            if sha(ROOT/name)==h:continue
            matching=[p for p in (DEST/'source_versions').glob('*.'+h+'.raw') if sha(p)==h]
            assert len(matching)==1,(name,h)
            differences.append(dict(path=name,sha256=h,archived=matching[0].relative_to(ROOT).as_posix()))
        tasks.append(dict(receipt=path.name,exitCode=d['exitCode'],lockedFiles=len(d['sourceSHA256']),changedSinceExecution=differences))
    preparations=[]
    for name in ['fixture','actual']:
        p=BASE/name/'mesh_process.json';d=read(p);assert d['status']=='terminal' and d['sourceStable'] and d['inputStable']
        for n,h in d['sourceHashes'].items():assert sha(Path(n))==h,n
        for n,h in d['outputSHA256'].items():assert sha(p.parent/n)==h,n
        preparations.append(dict(case=name,exitCode=d['exitCode'],seconds=d['elapsedSeconds']))
    upstream=[]
    with zipfile.ZipFile(str(ROOT/'tools/runtime/su2-thermal-build/su2.zip')) as z:
        for rel in ['Common/include/geometry/primal_grid/CPyramid.hpp','Common/src/geometry/primal_grid/CPyramid.cpp','Common/src/geometry/CPhysicalGeometry.cpp']:
            p=ROOT/'tools/runtime/su2-thermal-build/source'/rel
            assert p.read_bytes()==z.read(next(n for n in z.namelist() if n.endswith('/'+rel)))
            target=DEST/'upstream_source'/(p.name+'.'+sha(p)+'.raw');copy(p,target);upstream.append(dict(path=rel,sha256=sha(p)))
    covered={v['archived'] for v in copies}
    for p in sorted(DEST.rglob('*')):
        if p.is_file() and p.relative_to(ROOT).as_posix() not in covered:
            copies.append(dict(archived=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    manifest=dict(baseCommit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),files=copies,
        sourceBindings=sources,upstreamPyramidSources=upstream,matlabTasks=tasks,meshPreparations=preparations,
        nativeFlowRuns=1,actualRocketFlowRuns=0,newQualifiedArea=0,realCandidateAccepted=False,
        evidenceBytesBeforeManifest=sum(v['bytes'] for v in copies),nativeMouseOrScreenshotTest=False,
        analyticTestsReportedInSession=3,analyticStdoutReceiptAvailable=False,largeHistoricalResultsBackedUp=False,
        localOnlyCurrentFiles=[row for row in inventory if row['bytes']>=2*1024**2],
        unrelatedChangesNotCommitted=['web/viewer.js','web/section.js','tests/verifyWebSectionExport.py','docs/research/',
            'docs/validation/opencode_feature_retention/','docs/validation/opencode_section_export/','tests/testRocketFeatureRetentionMatrix.m','tests/testSectionExportReasons.js'])
    write(DEST/'delivery_manifest.json',manifest)
    print(json.dumps(dict(files=len(copies),evidenceBytes=manifest['evidenceBytesBeforeManifest'],outputFileLengthBytes=sum(unique.values()),matlabTasks=tasks),indent=2))
if __name__=='__main__':main()
