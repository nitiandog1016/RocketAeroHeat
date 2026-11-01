"""真实CGAL检查四项有独立几何真值的相交/非相交问题，不求流场。"""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from runCGALSurfaceAudit import audit_surface
def main():
    out=ROOT/'outputs/v136_validation/checker';assert not out.exists();out.mkdir(parents=True)
    cases=[('closed',[[0,0,0],[1,0,0],[0,1,0],[0,0,1]],[[0,2,1],[0,1,3],[1,2,3],[2,0,3]],0),
        ('shared_edge',[[0,0,0],[1,0,0],[0,1,0],[0,-1,0]],[[0,1,2],[1,0,3]],0),
        ('crossing',[[0,0,0],[2,0,0],[0,2,0],[.5,.5,-1],[.5,.5,1],[1.5,.5,0]],[[0,1,2],[3,4,5]],1),
        ('coplanar',[[0,0,0],[2,0,0],[0,2,0],[.2,.2,0],[.8,.2,0],[.2,.8,0]],[[0,1,2],[3,4,5]],1)]
    reports={}
    for name,p,f,n in cases:
        r=audit_surface(np.array(p,float),np.array(f,int),out/name,dict(kind='analytic-fixture',expectedPairs=n))
        assert r['pairCount']==n and r['originalCoordinatesAndOrientedFacesExactlyMatched'];reports[name]=r
    (out/'report.json').write_text(json.dumps(reports,indent=2)+'\n',encoding='utf8');print('4 actual CGAL fixtures passed')
if __name__=='__main__':main()
