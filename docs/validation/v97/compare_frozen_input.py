# -*- coding: utf-8 -*-
"""逐值比较原输入与小型冻结输入重算的全部MAT载荷；只排除容器头。"""
from pathlib import Path
import argparse
import json
import hashlib
import numpy as np
from scipy.io import loadmat


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('before',type=Path);parser.add_argument('after',type=Path)
    parser.add_argument('--report',type=Path,required=True);args=parser.parse_args();assert not args.report.exists()
    a=loadmat(args.before,simplify_cells=True);b=loadmat(args.after,simplify_cells=True);count=0

    def compare(x,y,path):
        nonlocal count
        if hasattr(x,'_fieldnames'):
            assert hasattr(y,'_fieldnames') and x._fieldnames==y._fieldnames,path
            for key in x._fieldnames:compare(getattr(x,key),getattr(y,key),path+'.'+key)
        elif isinstance(x,dict):
            assert isinstance(y,dict) and x.keys()==y.keys(),path
            for key in x:compare(x[key],y[key],path+'.'+key)
        elif isinstance(x,(list,tuple)) or isinstance(x,np.ndarray) and x.dtype==object:
            aa=list(x.flat) if isinstance(x,np.ndarray) else x
            bb=list(y.flat) if isinstance(y,np.ndarray) else y
            assert len(aa)==len(bb),path
            for i,(u,v) in enumerate(zip(aa,bb)):compare(u,v,path+'['+str(i)+']')
        else:
            try:same=np.array_equal(x,y,equal_nan=True)
            except TypeError:same=np.array_equal(x,y)
            assert same,path;count+=1

    compare({k:v for k,v in a.items() if not k.startswith('__')},
            {k:v for k,v in b.items() if not k.startswith('__')},'MAT')
    report=dict(status='passed-all-common-MAT-payload-objects-exactly-equal',savedObjectsExactlyEqual=count,
                excludedOnlyMATContainerMetadata=True,noNumericToleranceUsed=True,
                beforeSHA256=digest(args.before),afterSHA256=digest(args.after),scriptSHA256=digest(__file__))
    args.report.write_bytes((json.dumps(report,indent=2)+'\n').encode('utf8'));print(json.dumps(report))


if __name__=='__main__':main()
