"""模块职责：在内部原四面体状态上生成三维流线和九线横向度量；不授外缘或热资格。"""
import numpy as np
from scipy.spatial import cKDTree
from phaseAwareViscousAir import PhaseAwareViscousAir
from su2PhaseAirFields import inverse_energy


class TraceStop(ValueError):
    """可预期的物理/几何停止；调用者保存有效前缀和准确原因。"""
    def __init__(self,reason,point=None,face=None):
        super().__init__(reason);self.point=point;self.face=face


class TetraField:
    """模块一：原单元重心插值与面邻接行走；禁止把域外点拉回网格。"""
    def __init__(self, xyz, cells, conservative):
        self.x=np.asarray(xyz,dtype=float); self.cells=np.asarray(cells,dtype=np.int64)
        self.u=np.asarray(conservative,dtype=float); self.air=PhaseAwareViscousAir()
        if self.x.ndim!=2 or self.x.shape[1]!=3 or self.cells.ndim!=2 or self.cells.shape[1]!=4:
            raise ValueError('TUBE_MESH_SHAPE')
        if self.u.shape!=(len(self.x),5) or not np.isfinite(self.x).all() or not np.isfinite(self.u).all():
            raise ValueError('TUBE_STATE_SHAPE_OR_FINITE')
        if np.any(self.cells<0) or np.any(self.cells>=len(self.x)):raise ValueError('TUBE_NODE_RANGE')
        p=self.x[self.cells]; self.origin=p[:,0]
        mat=(p[:,1:]-p[:,0,None]).transpose(0,2,1)
        if np.any(np.linalg.det(mat)<=0):raise ValueError('TUBE_NONPOSITIVE_TET')
        self.inverse=np.linalg.inv(mat); self.lower=p.min(axis=1); self.upper=p.max(axis=1)
        self.tree=cKDTree(p.mean(axis=1));self.neighbor=np.full((len(p),4),-1,dtype=np.int64)
        half=np.concatenate([self.cells[:,[j for j in range(4) if j!=k]] for k in range(4)])
        owners=np.tile(np.arange(len(p)),4); opp=np.repeat(np.arange(4),len(p))
        keys,inv,count=np.unique(np.sort(half,axis=1),axis=0,return_inverse=True,return_counts=True)
        if np.any(count>2):raise ValueError('TUBE_NONMANIFOLD')
        order=np.argsort(inv,kind='stable');starts=np.r_[0,np.cumsum(count)[:-1]]
        internal=starts[count==2];a=order[internal];b=order[internal+1]
        self.neighbor[owners[a],opp[a]]=owners[b];self.neighbor[owners[b],opp[b]]=owners[a]
        # 仅用于包含关系的浮点容差，接受后不裁剪/归一化权重或移动坐标。
        self.bary_tol=2e-11; self.coord_tol=64*np.finfo(float).eps*max(1,float(abs(self.x).max()))
        self.states(self.u)  # 非法旧场直接拒绝，不能靠后处理修复。

    def weights(self,point,ids):
        v=np.einsum('nij,nj->ni',self.inverse[ids],point-self.origin[ids])
        return np.column_stack([1-v.sum(axis=1),v])

    def locate(self,point,hint=None):
        point=np.asarray(point,dtype=float)
        if point.shape!=(3,) or not np.isfinite(point).all():raise TraceStop('TUBE_POINT_INVALID')
        if hint is not None and 0<=hint<len(self.cells):
            seen=set();i=int(hint)
            for _ in range(64):
                if i in seen:break
                seen.add(i);w=self.weights(point,[i])[0]
                if w.min()>=-self.bary_tol:return i,w
                j=int(np.argmin(w));i=int(self.neighbor[i,j])
                if i<0:break
        # 最近质心只缩小首次查询；未包含则准确遍历包围盒，绝不最近邻补值。
        _,ids=self.tree.query(point,k=min(24,len(self.cells)));ids=np.atleast_1d(ids)
        w=self.weights(point,ids);good=np.flatnonzero(w.min(axis=1)>=-self.bary_tol)
        if len(good):j=good[0];return int(ids[j]),w[j]
        ids=np.flatnonzero(np.all(point>=self.lower-self.coord_tol,axis=1)&np.all(point<=self.upper+self.coord_tol,axis=1))
        if len(ids):
            w=self.weights(point,ids);good=np.flatnonzero(w.min(axis=1)>=-self.bary_tol)
            if len(good):j=good[0];return int(ids[j]),w[j]
        raise TraceStop('TUBE_OUTSIDE_ORIGINAL_DOMAIN')

    def sample(self,point,hint=None):
        i,w=self.locate(point,hint);u=w@self.u[self.cells[i]]
        return i,w,u

    def segment(self,start,end,hint=None):
        """连续线段沿实际面邻接行走，防止两端合法却跨过壁面或空洞。"""
        cell,_=self.locate(start,hint);at=0.;seen=set()
        for _ in range(128):
            a=self.weights(start,[cell])[0];b=self.weights(end,[cell])[0]
            if b.min()>=-self.bary_tol:return cell
            if cell in seen:raise TraceStop('TUBE_SEGMENT_AMBIGUOUS',start+at*(end-start))
            seen.add(cell);choices=[]
            for j in range(4):
                if b[j]<-self.bary_tol and a[j]-b[j]>0:
                    fraction=float(a[j]/(a[j]-b[j]))
                    if fraction>=at-2e-10:choices.append((fraction,j))
            if not choices:raise TraceStop('TUBE_SEGMENT_UNRESOLVED',start+at*(end-start))
            fraction,j=min(choices);fraction=max(at,fraction)
            point=start+fraction*(end-start);other=int(self.neighbor[cell,j])
            if other<0:raise TraceStop('TUBE_BOUNDARY_CROSSING',point,self.cells[cell,[k for k in range(4) if k!=j]])
            at=fraction;cell=other
        raise TraceStop('TUBE_SEGMENT_WALK_LIMIT',start+at*(end-start))

    def states(self,u):
        """模块二：同相态Air和原+300000能量坐标反算；不冻结H/S。"""
        u=np.asarray(u,dtype=float);rho=u[...,0]
        if np.any(rho<=0):raise TraceStop('TUBE_DENSITY_DOMAIN')
        v=u[...,1:4]/rho[...,None];speed=np.linalg.norm(v,axis=-1)
        e=u[...,4]/rho-.5*speed**2-300000
        try:
            temperature=inverse_energy(self.air,e);pressure=rho*self.air.R*temperature
            th=self.air.evaluate(pressure,temperature)
        except ValueError as exc:raise TraceStop('TUBE_THERMAL_DOMAIN: '+str(exc))
        return dict(rho=rho,velocity=v,speed=speed,p=pressure,T=temperature,
                    H=np.asarray(th['enthalpyJkg'])+.5*speed**2,S=np.asarray(th['entropyJkgK']))

    def direction(self,point,hint):
        i,w,u=self.sample(point,hint);state=self.states(u)
        if state['speed']<=1e-8:raise TraceStop('TUBE_STAGNATION')
        return state['velocity']/state['speed'],i


def trace(field,seed,ds,steps):
    """模块三：固定步长弧长RK4；每个实际阶段查域，失败时保留最后完整站。"""
    if not np.isfinite(ds) or ds<=0 or steps<1:raise ValueError('TUBE_STEP_CONTRACT')
    point=np.asarray(seed,dtype=float).copy();hint=None;points=[];cells=[];weights=[];states=[];reason='length-limit'
    failure_point=None;failure_face=None;phase='seed'
    for k in range(steps+1):
        try:
            hint,w,u=field.sample(point,hint);field.states(u)
            points.append(point.copy());cells.append(hint);weights.append(w);states.append(u)
            if k==steps:break
            a,_=field.direction(point,hint)
            phase='RK2-chord';trial=point+.5*ds*a;j=field.segment(point,trial,hint);b,j=field.direction(trial,j)
            phase='RK3-chord';trial=point+.5*ds*b;j=field.segment(point,trial,hint);c,j=field.direction(trial,j)
            phase='RK4-chord';trial=point+ds*c;j=field.segment(point,trial,hint);d,j=field.direction(trial,j)
            candidate=point+ds*(a+2*b+2*c+d)/6
            # 完整步端点也须合法；没有短步凑端点或夹回边界。
            phase='accepted-chord';j=field.segment(point,candidate,hint)
            _,_,end=field.sample(candidate,j);field.states(end);point=candidate
        except TraceStop as exc:
            reason=str(exc);failure_point=None if exc.point is None else np.asarray(exc.point).tolist()
            failure_face=None if exc.face is None else np.asarray(exc.face).tolist();break
    return dict(xyz=np.asarray(points).reshape(-1,3),cell=np.asarray(cells,dtype=np.int64),
                weights=np.asarray(weights).reshape(-1,4),U=np.asarray(states).reshape(-1,5),
                s=np.arange(len(points))*ds,stopReason=reason,requestedStations=steps+1,
                stopPhase=phase,stopPoint=failure_point,stopBoundaryNodes=failure_face)


LABELS=np.array([[0,0],[-1,0],[-.5,0],[.5,0],[1,0],[0,-1],[0,-.5],[0,.5],[0,1]],dtype=float)


def metrics(field,rays,delta):
    """模块四：两横向标签及流向组成真实三维管坐标；粗细差不冒充严格误差界。"""
    n=min(len(r['s']) for r in rays)
    if n<2:raise TraceStop('TUBE_NO_COMMON_PREFIX')
    p=np.asarray([r['xyz'][:n] for r in rays]);u=np.asarray([r['U'][:n] for r in rays]);s=field.states(u)
    tangent=s['velocity'][0]/s['speed'][0,:,None]
    fineB=(p[3]-p[2])/delta;coarseB=(p[4]-p[1])/(2*delta)
    fineE=(p[7]-p[6])/delta;coarseE=(p[8]-p[5])/(2*delta)
    def values(jb,je):
        b=jb-tangent*np.sum(jb*tangent,axis=1)[:,None]
        e=je-tangent*np.sum(je*tangent,axis=1)[:,None]
        area=np.einsum('ij,ij->i',tangent,np.cross(jb,je))
        return np.linalg.norm(b,axis=1),np.linalg.norm(e,axis=1),area
    hb,he,area=values(fineB,fineE);cb,ce,ca=values(coarseB,coarseE)
    valid=(area>1e-8)&(hb>1e-8)&(he>1e-8)
    # 聚焦后不允许通过后续恢复正号来重新授予标签一一映射。
    valid=np.logical_and.accumulate(valid)
    flux=s['rho'][0]*s['speed'][0]*area
    mass=flux/flux[0]
    h=s['H'][0];entropy=s['S'][0]
    return dict(commonStations=n,valid=valid,s=rays[0]['s'][:n],Jbeta=fineB,Jeta=fineE,
        hBeta=hb,hEta=he,signedArea=area,coarseHBeta=cb,coarseHEta=ce,coarseArea=ca,
        labelIndicator=np.maximum(np.linalg.norm(fineB-coarseB,axis=1),np.linalg.norm(fineE-coarseE,axis=1))/3,
        massFluxRatio=mass,HRelativeDrift=(h-h[0])/max(abs(h[0]),1.),entropyDrift=entropy-entropy[0],
        HBeta=(s['H'][3]-s['H'][2])/delta,HEta=(s['H'][7]-s['H'][6])/delta,
        SBeta=(s['S'][3]-s['S'][2])/delta,SEta=(s['S'][7]-s['S'][6])/delta,
        HAlong=np.gradient(h,rays[0]['s'][:n]),SAlong=np.gradient(entropy,rays[0]['s'][:n]))


def tube(field,seed,beta,eta,delta,ds,steps):
    """模块五：九条真实邻线一次生产；不替换失败标签、不提升研究场资格。"""
    seed=np.asarray(seed);beta=np.asarray(beta);eta=np.asarray(eta)
    if delta<=0 or not np.isfinite(delta):raise ValueError('TUBE_LABEL_SPACING')
    if beta.shape!=(3,) or eta.shape!=(3,) or abs(np.linalg.det(np.array([beta,eta,np.cross(beta,eta)])))<1e-8:
        raise ValueError('TUBE_LABEL_DIRECTIONS')
    rays=[trace(field,seed+delta*(a*beta+b*eta),ds,steps) for a,b in LABELS]
    return rays,metrics(field,rays,delta)
