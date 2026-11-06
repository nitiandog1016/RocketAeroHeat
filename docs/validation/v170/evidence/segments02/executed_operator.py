"""模块用途：原电极重叠测量段的面积积分与有标记的热参数归约；不伪造实验资格。"""
import csv
from pathlib import Path
import numpy as np

from naca3344Model import DATA,D,R,SHOULDER,END,source_data,heat_parameter,recovery_temperature
from naca3344Normalization import published_recovery
from su2ThermalRuntime import read,sha


def segments():
    """模块一：新段界绑定独立原图复核，同时保持旧20个观测值和横坐标。"""
    folder=DATA/'segments_v170';summary=read(folder/'summary.json')
    for name,h in summary['generated_csv_sha256'].items():assert sha(folder/name)==h
    for name,h in summary['old_inputs_sha256_before_and_after_identical'].items():
        if name.startswith('tests/'):
            assert sha(DATA.parents[2]/name)==h
    with (folder/'nominal_segments.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    assert len(rows)==20 and summary['nominal_segments_overlap']
    return rows,summary


def area(a,b):
    """模块二：名义球面/圆柱真实回转面积；此面积与离散三维壁面积分列。"""
    assert 0<=a<b<=END+1e-14
    a0=min(a,SHOULDER);b0=min(b,SHOULDER)
    # 和差化积避免鼻端两个接近1的cos相减。
    sphere=4*np.pi*R**2*np.sin((a0+b0)/(2*R))*np.sin((b0-a0)/(2*R))
    return sphere+2*np.pi*R*max(0,b-max(a,SHOULDER))


def integrate_profile(s,q_into_wall,provenance):
    """模块三：先积分有符号热流；输入须是有来源的周向平均，默认不推定轴对称。"""
    assert provenance.get('fieldSource') and provenance.get('circumferentialReduction')
    s=np.asarray(s,float);q=np.asarray(q_into_wall,float)
    assert s.ndim==q.ndim==1 and s.shape==q.shape and len(s)>=2
    assert np.isfinite(s).all() and np.isfinite(q).all() and np.all(np.diff(s)>0)
    assert s[0]>=0 and s[-1]<=END+1e-14
    rows,summary=segments();bounds=np.array([[float(r[k]) for k in ('s_left_m','s_right_m')] for r in rows])
    assert s[0]<=bounds.min() and s[-1]>=bounds.max(), '不能外插缺失的实测段热流'
    gx,gw=np.polynomial.legendre.leggauss(8)
    answer=[]
    for row,(a,b) in zip(rows,bounds):
        cuts=np.unique(np.r_[a,s[(s>a)&(s<b)],SHOULDER if a<SHOULDER<b else a,b])
        lo,hi=cuts[:-1],cuts[1:]
        x=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*gx/2
        radius=R*np.sin(np.minimum(x/R,np.pi/2))
        value=np.interp(x.ravel(),s,q).reshape(x.shape)
        power=float(np.sum(value*radius*2*np.pi*gw*(hi-lo)[:,None]/2))
        ar=float(area(a,b));assert abs(ar/float(row['area_full_circumference_m2'])-1)<2e-12
        answer.append(dict(row=int(row['heat_row_1based']),boundsM=[a,b],nominalFullAreaM2=ar,
            powerIntoWallW=power,meanQIntoWallWm2=power/ar,
            centerCandidatesS_m=dict(plot=float(row['original_plot_s_over_D'])*D,
                oppositeTemperatureStation=float(row['temperature_station_s_over_D'])*D,
                voltageIntervalMidpoint=float(row['voltage_interval_arc_midpoint_s_over_D'])*D),
            strictLocalAbscissaIdentityClosed=row['strict_local_abscissa_identity_closed']=='True',
            measuredParameter=float(row['original_Nu_s_over_sqrt_Re_s']),experimentalQualified=False))
    return dict(segments=answer,profileRepresentation='piecewise-linear signed circumferential-mean q(s)',
        provenance=provenance,nominalMeasuredUnionAreaM2=float(area(bounds.min(),bounds.max())),
        sumOverlappingAreasM2=float(sum(r['nominalFullAreaM2'] for r in answer)),
        overlappingPowersAreNotTotalPower=True,all20AbscissasClosed=summary['strict_all20_local_abscissa_identity_closed'],
        discreteWallPowerMapped=False,experimentalQualified=False)


def reduce_segment(segment,center_choice,given_tw,edge,edge_method):
    """模块四：按指定中心的同一局部外缘归约；平均热流不等于平均Nu/√Re。"""
    assert center_choice in segment['centerCandidatesS_m']
    assert edge_method and edge_method!='freestream-substituted-for-local-edge'
    c=source_data()[0]['freeStream'];s=segment['centerCandidatesS_m'][center_choice]
    cr=float(published_recovery(s));tr=float(recovery_temperature(edge['temperatureK'],c['totalTemperature_K'],cr,
        'NACA-TN3344-equation12-Figure3b-published-reduction-curve'))
    result=heat_parameter(s,segment['meanQIntoWallWm2'],given_tw,edge,tr,
        dict(edge=edge_method,recovery='published-curve; digitized'))
    result={k:float(v) if np.ndim(v)==0 and k!='experimentalQualification' else v for k,v in result.items()}
    result.update(row=segment['row'],centerChoice=center_choice,s_m=s,recoveryTemperatureK=tr,
        measuredParameter=segment['measuredParameter'],strictLocalAbscissaIdentityClosed=segment['strictLocalAbscissaIdentityClosed'],
        relativeDifference=float(result['parameter']/segment['measuredParameter']-1),
        originalTransportTableEquivalent=False,discreteWallPowerMapped=False,experimentalQualification=False,
        definition='finite-area mean q reduced with explicitly selected center state; no claim of full original reduction')
    return result


def transport_reduction_ratio(edge_a,edge_b,Tw,Tr_a,Tr_b):
    """模块五：固定同一s/q/Tw时两外缘模型的精确归约比；不把原表未知变成等价。"""
    for edge in (edge_a,edge_b):
        assert all(np.isfinite(edge[k]) and edge[k]>0 for k in ('densityKgM3','speedMS','viscosityPas','conductivityWmK'))
    assert np.isfinite(Tw+Tr_a+Tr_b) and Tw>Tr_a and Tw>Tr_b
    return float(np.sqrt(edge_b['viscosityPas']/edge_a['viscosityPas'])*edge_a['conductivityWmK']/edge_b['conductivityWmK']
        *np.sqrt(edge_a['densityKgM3']*edge_a['speedMS']/(edge_b['densityKgM3']*edge_b['speedMS']))*(Tw-Tr_a)/(Tw-Tr_b))
