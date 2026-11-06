"""原节点法向列的有限分辨诊断；不修改外缘判据、不补来流、不授实验资格。"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/'outputs/v171_validation/real01'
AUDIT = ROOT/'outputs/v171_validation/audit_real02'
RADIUS = .0508
DIAMETER = .1016
WINDOW = .001
STEPS = (1000, 1500)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024**2), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def variation_band(y, values):
    """模块一：以每区间绝对变化为权重，避免按节点累加梯度造成密网偏权。"""
    y = np.asarray(y, float)
    f = np.asarray(values, float)
    assert y.ndim == f.ndim == 1 and len(y) == len(f) and len(y) >= 2
    assert np.isfinite(y).all() and np.isfinite(f).all() and np.all(np.diff(y) > 0)
    dy, df = np.diff(y), np.diff(f)
    weight = abs(df)
    total = float(weight.sum())
    base = dict(nodeCount=len(y), intervalCount=len(y)-1, windowM=[float(y[0]), float(y[-1])],
                startValue=float(f[0]), endValue=float(f[-1]), minimumValue=float(f.min()),
                maximumValue=float(f.max()), totalAbsoluteVariation=total,
                signedChange=float(f[-1]-f[0]), observableVariation=bool(total > 0))
    if total == 0:
        return dict(base, band10to90=None, peakAbsoluteGradient=None)
    cdf = np.r_[0., np.cumsum(weight)]
    def crossing(fraction):
        target = fraction*total
        j = min(int(np.searchsorted(cdf, target, side='right'))-1, len(weight)-1)
        assert weight[j] > 0
        return float(y[j]+(target-cdf[j])/weight[j]*dy[j]), j
    low, lower = crossing(.1)
    high, upper = crossing(.9)
    j = int(np.argmax(weight/dy))
    interval = dict(leftIndex=j, rightIndex=j+1, boundsM=[float(y[j]), float(y[j+1])],
                    midpointM=float((y[j]+y[j+1])/2), widthM=float(dy[j]),
                    signedChange=float(df[j]), absoluteGradient=float(weight[j]/dy[j]),
                    totalVariationFraction=float(weight[j]/total))
    widths = dy[lower:upper+1]
    band = dict(boundsM=[low,high], widthM=high-low, leftIntervalIndex=lower,
                rightIntervalIndex=upper, intersectedIntervals=upper-lower+1,
                bracketingNodes=upper-lower+2,
                originalNodesInside=int(np.count_nonzero((y>=low)&(y<=high))),
                minimumIntervalM=float(widths.min()), maximumIntervalM=float(widths.max()),
                largestIntervalOverBandWidth=float(widths.max()/(high-low)))
    return dict(base, band10to90=band, peakAbsoluteGradient=interval)


def selected_stations():
    """模块二：仅固定原温度/电极锚点；后段不改变原实验中心身份。"""
    result = [('nose-axis', None, 0.)]
    for row in (1,4,8,11,14,15):
        result.append(('station-%02d'%row, row, RADIUS*np.deg2rad(6*row)))
    for row,offset in ((17,27),(20,67)):
        result.append(('station-%02d'%row,row,np.pi*RADIUS/2+offset/64*.0254))
    return result


def fixture(out):
    """模块三：不均匀网格的线性真值、往返变化与恒值检查，只核验诊断定义。"""
    start = time.monotonic()
    out.mkdir(parents=True, exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    y = np.array([0.,.01,.04,.15,.4,.7,1.])
    linear = variation_band(y, 17+8*y)
    assert abs(linear['totalAbsoluteVariation']-8) < 1e-13
    assert np.max(abs(np.array(linear['band10to90']['boundsM'])-[.1,.9])) < 1e-13
    folded = variation_band(np.array([0.,.5,1.]),np.array([0.,1.,0.]))
    assert folded['signedChange'] == 0 and folded['totalAbsoluteVariation'] == 2
    assert np.max(abs(np.array(folded['band10to90']['boundsM'])-[.1,.9])) < 1e-13
    constant = variation_band(y,np.full(len(y),3.))
    assert not constant['observableVariation'] and constant['band10to90'] is None
    rejected=[]
    for name,x,f in [('duplicate-distance',[0.,0.],[1.,2.]),('nonfinite',[0.,1.],[1.,np.nan])]:
        try: variation_band(x,f)
        except AssertionError: rejected.append(name)
        else: raise AssertionError('应拒绝：'+name)
    report=dict(passed=True,linearQuantileMaximumError=float(np.max(abs(np.array(linear['band10to90']['boundsM'])-[.1,.9]))),
                foldedTotalVariation=2.,constantHasNoDefinedBand=True,rejected=rejected,
                sourceSHA256=sha(__file__),elapsedSeconds=time.monotonic()-start,
                physicalQualification=False,experimentalQualification=False)
    save(out/'report.json',report)
    print(json.dumps(report,ensure_ascii=False))


def read_selected_state(path, expected_points):
    """只读需要的原生守恒量及T/p；不用最近邻、不重算CFD、不修改保存字段。"""
    with path.open(encoding='utf-8') as stream:
        names=next(csv.reader(stream))
    desired=['PointID','x','y','z','Density','Momentum_x','Momentum_y','Momentum_z','Pressure','Temperature']
    data=np.loadtxt(str(path),delimiter=',',skiprows=1,usecols=[names.index(k) for k in desired])
    assert data.shape==(len(expected_points),len(desired)) and np.isfinite(data).all()
    assert np.array_equal(data[:,0],np.arange(len(expected_points)))
    assert np.max(abs(data[:,1:4]-expected_points)) < 1e-14 and np.all(data[:,4]>0)
    return data[:,5:8]/data[:,4,None],data[:,8],data[:,9]


def scalar_ranges(values):
    return [float(min(values)),float(max(values))]


def review(out):
    """模块四：仅1000/1500两帧一次读取，全部原周向列与既存审核记录绑定。"""
    started=time.monotonic()
    out.mkdir(parents=True,exist_ok=False)
    (out/'executed_review.py').write_bytes(Path(__file__).read_bytes())
    audit=json.loads((AUDIT/'report.json').read_text(encoding='utf-8'))
    assert audit['savedFrameInterfacePassed'] and not audit['normalRunCompleted']
    paths=[Path(__file__).resolve(),ROOT/'scripts/validateIdealAir.py',AUDIT/'report.json',
           CASE/'geometry.json',CASE/'mesh_arrays.npz']
    for step in STEPS:
        paths.extend([CASE/('restart_flow_%06d.csv'%step),AUDIT/('measurement_%06d.json'%step)])
    locks={str(p):sha(p) for p in paths}
    for p in paths:
        if str(p) in audit['inputSHA256']:
            assert locks[str(p)]==audit['inputSHA256'][str(p)], '来源偏离已审核完整场'
    for step in STEPS:
        old=next(row for row in audit['frames'] if row['localIteration']==step)
        assert old['interfacePassed'] and locks[str(AUDIT/old['measurementJSON'])]==old['measurementSHA256']
    save(out/'input_manifest.json',locks)
    meta=json.loads((CASE/'geometry.json').read_text(encoding='utf-8'))
    with np.load(str(CASE/'mesh_arrays.npz')) as arrays:
        points=arrays['points'].copy()
        ix=arrays['columnNodeIDs'].copy()
    y=np.asarray(meta['offsets']);s=np.asarray(meta['s']);az=np.asarray(meta['azimuth'])
    assert ix.shape==(len(y),len(s),len(az)) and np.all(np.diff(y)>0)
    cut=int(np.searchsorted(y,WINDOW,side='right')-1)
    assert 2<=cut<len(y)-2
    selection=[]
    for name,row,requested in selected_stations():
        matches=np.flatnonzero(abs(s-requested)<1e-13)
        assert len(matches)==1,'必须是唯一原列站，不允许近邻替代'
        i=int(matches[0]);mapping={};columns=[]
        for k,phi in enumerate(az):
            ids=ix[:,i,k]
            key=tuple(int(v) for v in ids)
            if key not in mapping:
                mapping[key]=len(columns)
                theta=min(s[i]/RADIUS,np.pi/2)
                normal=np.array([np.cos(theta),np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi)])
                displacement=points[ids]-points[ids[0]]
                assert np.max(abs(displacement-y[:,None]*normal))<1e-13
                columns.append(dict(uniqueColumn=len(columns),originalNodeIDs=ids.tolist(),azimuthIndices=[],normal=normal.tolist()))
            columns[mapping[key]]['azimuthIndices'].append(k)
        selection.append(dict(name=name,originalTemperatureStationRow=row,sIndex=i,s_m=float(s[i]),s_over_D=float(s[i]/DIAMETER),
                              requestedM=requested,azimuthSlots=len(az),independentColumns=len(columns),columns=columns))
    all_frames={}; summaries=[]; profile_rows=[]
    for step in STEPS:
        velocity,pressure,temperature=read_selected_state(CASE/('restart_flow_%06d.csv'%step),points)
        old=json.loads((AUDIT/('measurement_%06d.json'%step)).read_text(encoding='utf-8'))
        old_segments={int(row['row']):row for row in old['segments']}
        frame=[];station_summaries=[]
        for station in selection:
            station_results=[]
            for column in station['columns']:
                ids=np.asarray(column['originalNodeIDs']);normal=np.asarray(column['normal'])
                local=velocity[ids];un=local@normal
                tangent=local-un[:,None]*normal
                ut=np.linalg.norm(tangent,axis=1)
                fields=dict(temperatureK=temperature[ids],pressurePa=pressure[ids],normalVelocityMS=un,
                            tangentialSpeedMS=ut,speedMS=np.linalg.norm(local,axis=1))
                prior=[]
                for k in column['azimuthIndices']:
                    if station['originalTemperatureStationRow'] is None:
                        prior.append(dict(azimuthIndex=k,status='stagnation-limit-not-defined',accepted=False))
                    else:
                        entry=old_segments[station['originalTemperatureStationRow']]['centerColumns']['oppositeTemperatureStation'][k]
                        assert entry['originalNodeIDs']==ids.tolist() and entry['sIndex']==station['sIndex']
                        prior.append(dict(azimuthIndex=k,status=entry['status'],accepted=entry['accepted']))
                metrics={}
                for name,f in fields.items():
                    metrics[name]=dict(fullColumn=variation_band(y,f),nearWallFixedWindow=variation_band(y[:cut+1],f[:cut+1]),
                                       outerFixedWindow=variation_band(y[cut:],f[cut:]))
                    total=metrics[name]['fullColumn']['totalAbsoluteVariation']
                    metrics[name]['wallFirstIntervalVariationFraction']=float(abs(f[1]-f[0])/total) if total>0 else None
                row=dict(station=station['name'],sIndex=station['sIndex'],uniqueColumn=column['uniqueColumn'],
                         azimuthIndices=column['azimuthIndices'],normal=column['normal'],originalNodeIDs=ids.tolist(),
                         originalEdgeResults=prior,metrics=metrics,physicalQualification=False)
                station_results.append(row);frame.append(row)
                for j,node in enumerate(ids):
                    profile_rows.append([step,station['name'],station['sIndex'],column['uniqueColumn'],
                                         '|'.join(map(str,column['azimuthIndices'])),j,int(node),float(y[j]),
                                         *[float(fields[key][j]) for key in ('temperatureK','pressurePa','normalVelocityMS','tangentialSpeedMS','speedMS')]])
            summary=dict(name=station['name'],s_over_D=station['s_over_D'],independentColumns=len(station_results),
                         originalAcceptedEdges=sum(int(e['accepted']) for r in station_results for e in r['originalEdgeResults']))
            for field,window,prefix in [('temperatureK','nearWallFixedWindow','nearWallTemperature'),('pressurePa','outerFixedWindow','outerPressure'),
                                        ('temperatureK','outerFixedWindow','outerTemperature')]:
                metrics=[r['metrics'][field][window] for r in station_results]
                bands=[m['band10to90'] for m in metrics if m['band10to90'] is not None]
                peaks=[m['peakAbsoluteGradient'] for m in metrics if m['peakAbsoluteGradient'] is not None]
                summary[prefix]=dict(definedColumns=len(bands),bandWidthRangeM=scalar_ranges([b['widthM'] for b in bands]) if bands else None,
                    bracketingNodesRange=scalar_ranges([b['bracketingNodes'] for b in bands]) if bands else None,
                    maximumIntervalRangeM=scalar_ranges([b['maximumIntervalM'] for b in bands]) if bands else None,
                    peakMidpointRangeM=scalar_ranges([p['midpointM'] for p in peaks]) if peaks else None,
                    peakIntervalWidthRangeM=scalar_ranges([p['widthM'] for p in peaks]) if peaks else None,
                    totalVariationRange=scalar_ranges([m['totalAbsoluteVariation'] for m in metrics]))
            summary['temperatureFirstIntervalFractionRange']=scalar_ranges([r['metrics']['temperatureK']['wallFirstIntervalVariationFraction'] for r in station_results])
            station_summaries.append(summary)
        all_frames[step]=frame
        save(out/('columns_%06d.json'%step),frame)
        summaries.append(dict(localIteration=step,stations=station_summaries))
    movements=[]
    for a,b in zip(all_frames[1000],all_frames[1500]):
        assert a['station']==b['station'] and a['originalNodeIDs']==b['originalNodeIDs']
        changes={}
        for field in a['metrics']:
            changes[field]={}
            for window in ('fullColumn','nearWallFixedWindow','outerFixedWindow'):
                ma,mb=a['metrics'][field][window],b['metrics'][field][window]
                ba,bb=ma['band10to90'],mb['band10to90']
                pa,pb=ma['peakAbsoluteGradient'],mb['peakAbsoluteGradient']
                changes[field][window]=dict(bandBoundsMovementM=[float(x-y0) for x,y0 in zip(bb['boundsM'],ba['boundsM'])] if ba and bb else None,
                    peakMidpointMovementM=pb['midpointM']-pa['midpointM'] if pa and pb else None,
                    peakIntervalIndexChange=pb['leftIndex']-pa['leftIndex'] if pa and pb else None,
                    totalVariationChange=mb['totalAbsoluteVariation']-ma['totalAbsoluteVariation'])
        movements.append(dict(station=a['station'],uniqueColumn=a['uniqueColumn'],azimuthIndices=a['azimuthIndices'],changes=changes))
    save(out/'cross_frame_movements.json',movements)
    with (out/'selected_original_profiles.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['localIteration','station','sIndex','uniqueColumn','azimuthIndices','normalIndex','PointID','distanceM',
                         'temperatureK','pressurePa','normalVelocityMS','tangentialSpeedMS','speedMS'])
        writer.writerows(profile_rows)
    stable={name:sha(Path(name))==h for name,h in locks.items()}
    assert all(stable.values())
    report=dict(diagnosticCompleted=True,baselineCommit='2865300',savedIterations=list(STEPS),
        sourceTaskNormalCompleted=False,sourceSavedFrameAuditPassed=True,sourceTaskTimedOut=True,
        stationCount=len(selection),azimuthSlotsPerStation=len(az),normalNodesPerColumn=len(y),
        totalAzimuthColumnSlots=sum(r['azimuthSlots'] for r in selection),independentColumnsPerFrame=sum(r['independentColumns'] for r in selection),
        nearWallWindowRequestedM=WINDOW,nearWallWindowOriginalEndIndex=cut,nearWallWindowActualEndM=float(y[cut]),
        outerWindowStartM=float(y[cut]),selection=selection,frames=summaries,
        bandDefinition='10-90 percent of sum(abs(delta f)); piecewise-linear interval variation, never node-sum gradient',
        windowMeaning='fixed geometric diagnostic windows; 1mm follows documented initial blending scale and is not a boundary-layer edge',
        outerPressureMeaning='observable pressure-transition candidate, not an automatically identified bow shock',
        originalEdgePolicyChanged=False,freestreamFallbackUsed=False,newDelta99Assigned=False,
        physicalQualification=False,gridAndDomainQualified=False,experimentalQualification=False,
        inputSHA256=locks,inputsStable=stable,elapsedSeconds=time.monotonic()-started)
    save(out/'report.json',report)
    print(json.dumps({k:v for k,v in report.items() if k not in {'selection','frames','inputSHA256','inputsStable'}},ensure_ascii=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--fixture',action='store_true')
    args=parser.parse_args()
    (fixture if args.fixture else review)(Path(args.output).resolve())
