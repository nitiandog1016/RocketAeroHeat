"""模块职责：保留首次越壁轨迹，量化守卫恢复的影响；不重积分或重算流场。"""
import numpy as np
from runSU2RocketPhaseV222 import ROOT,read,save,sha
from phaseEulerRegionalVolume import read_tetra_mesh
from auditForebodyStreamtubesV241 import segment_boundary_hits


def compare():
    base=ROOT/'outputs/v241_validation';old=read(base/'streamtubes01/manifest.json');new=read(base/'streamtubes02/manifest.json')
    x,_,markers=read_tetra_mesh(ROOT/'outputs/v240_validation/run01/mesh.su2');rows=[];identical=0
    for before,after in zip(old['tubes'],new['tubes']):
        a=np.load(str(base/'streamtubes01'/before['file']));b=np.load(str(base/'streamtubes02'/after['file']))
        n=after['commonStations'];assert n==before['commonStations']
        for j in range(9):
            for key in ('xyz','cell','weights','U','s'):
                assert np.array_equal(a['ray%d_%s'%(j,key)][:n],b['ray%d_%s'%(j,key)][:n])
            identical+=n
        for key in a.files:
            if key.startswith('metric_'):assert np.array_equal(a[key],b[key])
        crossings=segment_boundary_hits(a['ray0_xyz'],x,markers)
        rows.append(dict(tube=after['id'],originalCenterStations=len(a['ray0_s']),correctedCenterStations=len(b['ray0_s']),
           originalCenterInteriorBoundaryCrossings=crossings,commonNineRayPrefixBitEqual=True,allSavedMetricsBitEqual=True))
    assert len(rows[1]['originalCenterInteriorBoundaryCrossings'])==1
    assert rows[1]['originalCenterStations']==274 and rows[1]['correctedCenterStations']==135
    save(base/'recovery01/report.json',dict(passed=True,tubes=rows,bitEqualCommonRayStations=identical,
      firstAuditFailureAttributedToTube=2,firstAuditFailureSegmentZeroBased=134,
      firstCommentaryDidNotIdentifyTube=True,planInitialTubeNumberCorrection='first->second; original plan snapshot retained',
      wholeOriginalRawTrajectoriesValid=False,allNewSavedSegmentsPreviouslyAudited=True,
      originalCFDFieldModified=False,flowUpdates=0,trajectoryReintegrations=0,
      unmeteredTargetedReadOnlyFollowups=3,physicalQualification=False))


if __name__=='__main__':compare()
