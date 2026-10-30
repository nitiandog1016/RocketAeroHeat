# -*- coding: utf-8 -*-
"""独立复核七参数实际拟合及未进入拟合的36个检查点，不把拟合通过当成区域通过。"""
from pathlib import Path
import argparse
import json

import numpy as np
from scipy.io import loadmat

from auditEllipticConicShockLayer import reference, geometric_check
from auditPressureDrivenSurfaceFlow import rows
from auditThermalPressureSurfaceFlow import sha

ROOT=Path(__file__).resolve().parents[1]


def audit(evidence,process,target):
    assert not target.exists();s=loadmat(evidence,simplify_cells=True)
    p=json.loads(process.read_text(encoding='utf-8-sig'))
    assert p['status']=='terminal' and p['exitCode']==0 and all(p['sourceStableAfter'].values())
    fs=s['fs'];up={'pressure':fs['pressure'],'temperature':fs['temperature'],
                 'density':fs['density'],'velocityVector':fs['flowVelocityVector']}
    fit=s['bodyFit'];checks=[];residuals={}
    for key,result in [('initial',fit['initial']),('result',fit['result']),('heldOut',s['heldOut'])]:
        result['columns']=rows(result['columns'])
        for c in result['columns']:checks.append(reference(c,up,c['shock']['gas']))
        residuals[key]=geometric_check(s['geometry'],result)
    assert np.linalg.norm(residuals['result'])<np.linalg.norm(residuals['initial'])
    assert len(residuals['result'])==9 and len(residuals['heldOut'])==36
    assert not s['heldOut']['usedDuringFit']
    for name in ['initial','result','heldOut']:
        result=s['heldOut'] if name=='heldOut' else fit[name]
        np.testing.assert_allclose(result['residualM'],residuals[name],rtol=0,atol=1e-11)
    tolerance=fit['bodyMatchToleranceM'];matched=np.max(abs(residuals['result']))<=tolerance
    held=np.max(abs(residuals['heldOut']))<=tolerance
    assert bool(fit['shapeMatchingConverged'])==bool(matched)
    assert bool(s['report']['heldOutBodyMatchingPassed'])==bool(held)
    assert not fit['boundaryLayerEdgeAvailable'] and not fit['productionQualified']
    sources=0
    for rel,digest in p['sourceSHA256'].items():
        if rel.startswith('outputs/'):continue
        file=ROOT/rel;assert file.exists() and sha(file)==digest,rel;sources+=1
    report=dict(status='passed-independent-all-fitted-and-held-out-physical-column-audit',
                evidenceSHA256=sha(evidence),processSHA256=sha(process),columnsChecked=len(checks),
                fittedColumns=9,heldOutColumns=36,
                maximumRelativeStateDifference=max(c['maximumRelativeStateDifference'] for c in checks),
                maximumNormalDistanceDifferenceM=max(c['maximumDistanceDifferenceM'] for c in checks),
                maximumFitResidualM=float(np.max(abs(residuals['result']))),
                maximumHeldOutResidualM=float(np.max(abs(residuals['heldOut']))),
                bodyMatchToleranceM=float(tolerance),sampledBodyMatchingConverged=bool(matched),
                heldOutBodyMatchingPassed=bool(held),sourceSHARecordsVerified=sources,
                stagnationGlobalMatchingSolved=False,boundaryLayerEdgeAvailable=False,
                experimentalValidation=False,newProductionSolvedAreaM2=0)
    target.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('process',type=Path)
    p.add_argument('--report',required=True,type=Path);a=p.parse_args();audit(a.evidence,a.process,a.report)
