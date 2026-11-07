"""固定v187完整域的资源下界与启动拒绝；不把下界当作运行峰值预测。"""
from pathlib import Path
import argparse
import hashlib
import json
import time

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_SHA = 'd904a4079838eefcfe9f7e424b6dd2bba5c88e671b19ca796c8fc645622ae08b'
MESH_SHA = '00f96f3a0cb87a11c4218796cd1f5f25aebffb453c794e0cf97f7bd4abde3181'
SOURCE_COMMIT = 'f53f4609a0331d48efee69a4e4f3c3598378cc03'
COUNTS = dict(cells=705490, ghosts=26324, faces=1741600, vertices=335726,
              tetrahedra=70574, prisms=634916)
SIZES = dict(FlowGradients=200, LSQInterpGradients=824, WLSQGradWorkspace=312,
             LSQInterpWorkspace=2688, FlowState=288, FluidFVCell=680,
             FVInterface=512, Matrix=48, Vector3=24, Slice=16, Scalar=8)
BINS = (16, 32, 48, 64, 96, 128, 176, 256, 368, 512, 816, 1024, 1360, 2048)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def lower_bound():
    """逐项只计独立且初始化完成时并存的分配；共享视图只计切片描述符。"""
    n, g, f, v = (COUNTS[k] for k in ('cells', 'ghosts', 'faces', 'vertices'))
    c = n + g
    cv = 4 * COUNTS['tetrahedra'] + 6 * COUNTS['prisms']
    cf = 4 * COUNTS['tetrahedra'] + 5 * COUNTS['prisms']
    fv = 3 * f  # 每面至少三顶点；实际四边形只取三点，保持严格下界。
    rows = []

    def add(name, value, formula, source, group='payload'):
        rows.append(dict(name=name, bytes=value, formula=formula, source=source, group=group))

    # 模块一：展开v187已确认的必要载荷，避免将旧合计当作不可检查的常数。
    for name, unit, formula in [
        ('stageArrays', 8*7*8, 'C*8*7*8'), ('sourceTerms', 7*8, 'C*7*8'),
        ('cellFlowStates', 288, 'C*288'), ('cellGradients', 200+24, 'C*(200+24)'),
        ('cellGradientWorkspace', 312, 'C*312'), ('cellInterpolationWorkspace', 2688, 'C*2688'),
        ('cellInterpolationGradients', 824+48, 'C*(824+48)'), ('cellObjects', 680, 'C*680'),
        ('cellStageViewsAndUdf', 2*4*16+56, 'C*(2*4*16+56)'),
        ('cellGeometry', 24+8+8+24, 'C*(24+8+8+24)')]:
        add(name, c*unit, formula, 'fluidblock/ufluidblock/fluidfvcell + v187 ABI', 'v187')
    add('faceObjects', f*512, 'F*512', 'fvinterface + v187 ABI', 'v187')
    add('faceFlowAndGradients', f*(288+200+24+312), 'F*(288+200+24+312)', 'fluidblock:333-338', 'v187')
    add('faceCorePriorLowerBound', f*226, 'F*226', 'v187 complete_capacity.py; dL/dR each counted as 8', 'v187')
    # 模块二：独立拓扑、顶点与邻域。类的字段和采用未填充下界，随后按GC规格取整。
    add('faceVectorRemainder', f*32, 'F*(2*24-2*8)', 'coredata:104-105')
    add('faceStencil', f*48, 'F*(6*8)', 'coredata:99,113; fluidblock:326')
    add('faceInterpolation', f*72, 'F*(9*8)', 'onedinterp:78-87; fluidblock:330')
    add('vertices', v*(137+48+200+24+8), 'V*(137+2*24+200+24+8)', 'fvvertex:26-55; ufluidblock:307-310')
    add('gridObjects', f*65+n*76, 'F*65+N*76', 'usgrid:154-165,218-224; includes Object header 16')
    add('gridConnectivity', fv*8+cv*8+cf*12, 'FV*8+CV*8+CF*(8+4)', 'usgrid:167-193,245-269')
    add('gridTopArrays', f*8+n*8+v*24, 'F*8+N*8+V*24', 'usgrid:424-432; grid.vertices')
    add('cellTopology', cv*16+cf*24+cf*24, 'CV*(8+8)+CF*(8+8+4+4)+CF*24', 'ufluidblock:364-405')
    add('faceVertices', fv*8, 'FV*8', 'ufluidblock:357-360')
    add('vertexIncidence', (cv+fv)*8+v*32, '(CV+FV)*8+V*2*16', 'ufluidblock:409-415')
    add('cellOuterSlices', n*16*5, 'N*5*16', 'ufluidblock:328-331,362-363; coredata:24-38')
    add('interpolationCloud', (n+cf)*8+cf*8, '(N+CF)*8+CF*8', 'ufluidblock:588-606; all ghosts contain data')
    add('viscousCloud', (n+cf)*16, '(N+CF)*2*8', 'ufluidblock:623-634; pointers only')
    add('gasDynamic', (c+f)*16, '(C+F)*2*8', 'gas_state:52-65')
    add('faceDiffusion', f*32, 'F*4*8', 'fvinterface:136-140; multi_species_gas retained')
    add('denseCellMisc', n*34+c*32, 'N*(8+8+8+8+1+1)+C*2*16', 'fluidblock:258-267; coredata:22-38')
    add('solverObjectArrays', n*8+f*8, '(N+F)*8', 'ufluidblock:341-354; ghost refs excluded')
    # 模块三：只加入已知的小分配最小取整差，不计数组元数据、扩容、页空隙或GC预留。
    def rounding(name, count, payload, source):
        minimum = next(b for b in BINS if b >= payload)
        add(name, count*(minimum-payload), '%d*(%d-%d)' % (count, minimum, payload), source, 'rounding')
    rounding('cellClassRounding', c, 680, 'FluidFVCell + GC bins')
    rounding('vertexClassRounding', v, 137, 'FVVertex unpadded fields + GC bins')
    rounding('gridFaceClassRounding', f, 65, 'USGFace unpadded fields + GC bins')
    rounding('gridCellClassRounding', n, 76, 'USGCell unpadded fields + GC bins')
    rounding('gasArrayRounding', 2*(c+f), 8, 'massf/rho_s separately allocated')
    rounding('diffusionArrayRounding', 4*f, 8, 'jx/jy/jz/hs separately allocated')
    rounding('interpVectorRounding', c, 24, 'LSQInterpGradients.rho_s')
    rounding('interpScalarRounding', 3*c, 8, 'LSQInterpGradients.rho_sPhi/Max/Min')
    rounding('gradientArrayRounding', c+f+v, 24, 'FlowGradients.massf')
    rounding('udfArrayRounding', c, 56, 'FluidFVCell.Qudf')
    rounding('vertexGradientRounding', v, 200, 'FVVertex new FlowGradients')
    return dict(rows=rows, counts=dict(C=c, N=n, F=f, V=v, CV=cv, CF=cf, FV=fv),
                payloadLowerBoundBytes=sum(r['bytes'] for r in rows if r['group'] != 'rounding'),
                minimumAllocationLowerBoundBytes=sum(r['bytes'] for r in rows),
                fullRunPeakBytes=None, fitProven=False, physicalQualification=False)


def assess_budget(mesh, package_sha, memory_mib):
    """仅对准确包和原网格作必要条件拒绝；更大资源不等于允许启动。"""
    if package_sha != PACKAGE_SHA or mesh.get('sourceSHA256') != MESH_SHA:
        raise ValueError('资源下界只适用于已绑定的v187运行包和原完整域')
    required = dict(cells=705490, nodes=335726, wallFaces=19752, farfieldFaces=6572,
                    cellTypeCounts={'10':70574, '12':0, '13':634916, '14':0})
    if any(mesh.get(k) != value for k, value in required.items()):
        raise ValueError('资源核算网格计数或单元类型不符')
    if not 1.9646292565713 < mesh.get('wallAreaM2', 0) < 1.9646292565718:
        raise ValueError('资源核算原壁面积不符')
    if type(memory_mib) is not int or memory_mib <= 0:
        raise ValueError('资源预算须为正整数MiB')
    bound = lower_bound()['minimumAllocationLowerBoundBytes']
    limit = memory_mib * 1024**2
    return dict(status='impossible' if bound >= limit else 'not-excluded-unproven',
                lowerBoundBytes=bound, budgetBytes=limit, runAuthorized=False,
                reason='初始化必要分配下界；尚未包含完整运行峰值')


def inspect(out):
    """绑定已交付ABI、准确覆盖源和原配置；只保存小证据，不读取体场或启动引擎。"""
    start = time.monotonic()
    out = Path(out); out.mkdir(parents=True, exist_ok=False)
    result = dict(passed=False, nativeRuns=0, builds=0)
    try:
        base = ROOT/'outputs/v187_validation'
        paths = dict(scope=base/'scope03/report.json', prior=base/'capacity02/report.json',
                     config=base/'full/real01/lmrsim/config', manifest=base/'full/real01/case_manifest.json',
                     runtime=ROOT/'tools/runtime/eilmer-explicit-workspace-v1/runtime.json')
        docs = {k:json.loads(p.read_text()) for k,p in paths.items()}
        assert sha(paths['runtime']) == PACKAGE_SHA
        scope, cfg, runtime = docs['scope'], docs['config'], docs['runtime']
        assert scope['passed'] and scope['capacity']['sizes'] == SIZES
        assert runtime['sourceCommit'] == SOURCE_COMMIT
        assert all(flag in scope['compileCommand'] for flag in ['-d-version=newton_krylov', '-d-version=multi_species_gas'])
        assert cfg['interpolation_order'] == 2 and not cfg['use_extended_stencil'] and cfg['viscous']
        assert cfg['spatial_deriv_locn'] == 'cells' and cfg['block_0']['nfaces'] == COUNTS['faces']
        assert all(cfg['block_0']['boundary_'+str(i)]['ghost_cell_data_available'] for i in range(2))
        result.update(lower_bound())
        assert sum(r['bytes'] for r in result['rows'] if r['group']=='v187') == docs['prior']['partialNecessaryLowerBoundBytes']
        result['resourceDecision'] = assess_budget(docs['manifest']['mesh'], PACKAGE_SHA, 8192)
        src = Path('/opt/rocket-aeroheat/eilmer-v5.0.0/source')
        overlay = ROOT/'tools/runtime/eilmer-explicit-workspace-v1/overlay'
        specs = {
            'src/lmr/fluidblock.d':[(237,355)], 'src/lmr/ufluidblock.d':[(283,665)],
            'src/lmr/coredata.d':[(21,53),(98,121)], 'src/lmr/onedinterp.d':[(78,101)],
            'src/lmr/fvvertex.d':[(26,59)], 'src/geom/grid/usgrid.d':[(154,193),(218,270),(424,432)],
            'src/lmr/fvinterface.d':[(112,141)], 'src/lmr/fluidfvcell.d':[(179,253)],
            'src/gas/gas_state.d':[(52,65)], 'src/lmr/flowgradients.d':[(31,110)],
            'src/lmr/lsqinterp.d':[(23,34),(250,315)],
            'src/eilmer/ufluidblock.d':[(296,330)], 'src/eilmer/lsqinterp.d':[(23,34)]}
        hashes, excerpts = {}, []
        for name, ranges in specs.items():
            p = overlay/name if (overlay/name).exists() else src/name
            digest = sha(p)
            expected = runtime['filesSHA256'].get(str(p), scope['originalSHA256'].get(str(p)))
            if expected is not None:
                assert digest == expected, str(p)
            elif not name.startswith('src/eilmer/'):
                raise AssertionError('source not bound to build: '+name)
            hashes[str(p)] = digest
            lines = p.read_text().splitlines()
            excerpts.append(dict(path=str(p), sha256=digest, ranges=[dict(start=a,end=b,
                text='\n'.join('%d: %s' % (i,lines[i-1]) for i in range(a,b+1))) for a,b in ranges]))
            if p == overlay/name:
                original = src/name
                hashes[str(original)] = sha(original)
                assert hashes[str(original)] == scope['originalSHA256'][str(original)]
                original_lines = original.read_text().splitlines()
                excerpts.append(dict(path=str(original),sha256=sha(original),ranges=[dict(start=a,end=b,
                    text='\n'.join('%d: %s' % (i,original_lines[i-1]) for i in range(a,b+1))) for a,b in ranges]))
        gc = src.parent/'ldc2-1.41.0-linux-x86_64/import/core/internal/gc/impl/conservative/gc.d'
        hashes[str(gc)] = sha(gc); lines = gc.read_text().splitlines()
        assert '[ 16, 32, 48, 64, 96, 128, 176, 256, 368, 512, 816, 1024, 1360, 2048, 4096 ]' in lines[1681]
        excerpts.append(dict(path=str(gc),sha256=sha(gc),ranges=[dict(start=a,end=b,
            text='\n'.join('%d: %s' % (i,lines[i-1]) for i in range(a,b+1))) for a,b in [(1646,1683),(2158,2212)]]))
        hashes.update({str(p):sha(p) for p in paths.values()}); hashes[str(Path(__file__).resolve())] = sha(__file__)
        result.update(passed=True,inputSHA256=hashes,sourceCommit=SOURCE_COMMIT,
                      excludes=['array metadata/growth', 'GC pool slack and page gaps', 'temporary/IO buffers',
                                'runtime/mappings/stacks', 'later solver allocations'],
                      interpretation='Necessary simultaneously live allocations at completed initialization; not failed RSS, peak virtual size, or exact OOM site.')
        (out/'source_excerpts.json').write_text(json.dumps(excerpts,indent=2)+'\n')
        result['inputsStable'] = all(sha(p)==h for p,h in hashes.items())
        assert result['inputsStable']
    except Exception as exc:
        result['error'] = repr(exc)
        raise
    finally:
        result['elapsedSeconds'] = time.monotonic()-start
        (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('passed','minimumAllocationLowerBoundBytes','elapsedSeconds')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out',required=True)
    inspect(parser.parse_args().out)
