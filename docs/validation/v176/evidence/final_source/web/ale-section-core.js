/* 实际 X 截面与三角面的交线；沿交线保留面元常值及未解状态。 */
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.RocketSection=api;})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  const tau=2*Math.PI,finite=Number.isFinite;
  const distance=(a,b)=>Math.hypot(...a.map((v,j)=>v-b[j]));
  const angle=(p,c)=>((Math.atan2(p[2]-c[2],p[1]-c[1])+tau)%tau)*180/Math.PI;
  function geometry(vertices,faces,x){
    if(!finite(x))throw Error('截面 X 必须为有限数值。');
    let lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity],magnitude=0;
    vertices.forEach(v=>v.forEach((q,j)=>{lo[j]=Math.min(lo[j],q);hi[j]=Math.max(hi[j],q);magnitude=Math.max(magnitude,Math.abs(q));}));
    const size=Math.max(...hi.map((q,j)=>q-lo[j])),eps=Math.max(size*1e-10,magnitude*Number.EPSILON*128);
    if(!(size>0)||x<lo[0]-eps||x>hi[0]+eps)return {x,segments:[],components:[],coplanarFaces:0};
    // 截面与网格环重合时采用 +X 侧面元常值；最前端采用 -X 侧。
    const trace=x>=hi[0]-eps?-1:1,segments=[];let coplanarFaces=0;
    faces.forEach((face,index)=>{
      const vs=face.map(i=>vertices[i]),ds=vs.map(v=>Math.abs(v[0]-x)<=eps?0:v[0]-x);
      if(ds.every(d=>d===0)){coplanarFaces++;return;}
      if(ds.every(d=>d>0)||ds.every(d=>d<0))return;
      if(ds.filter(d=>d===0).length===2&&ds.find(d=>d!==0)*trace<0)return;
      const points=[],pointKeys=[];
      function add(p,key){if(!points.some(q=>distance(p,q)<=eps)){points.push(p);pointKeys.push(key);}}
      for(let j=0;j<3;j++){
        const k=(j+1)%3;if(ds[j]===0)add([x,vs[j][1],vs[j][2]],'v:'+face[j]);
        if(ds[j]*ds[k]<0){const t=(x-vs[j][0])/(vs[k][0]-vs[j][0]);add([x,vs[j][1]+t*(vs[k][1]-vs[j][1]),vs[j][2]+t*(vs[k][2]-vs[j][2])],'e:'+Math.min(face[j],face[k])+':'+Math.max(face[j],face[k]));}
      }
      if(points.length!==2)return;
      const length=distance(points[0],points[1]);if(length<=eps)return;
      segments.push({face:index,a:points[0],b:points[1],length,_keys:pointKeys});
    });
    // 模块二：共享顶点/边采用拓扑键，避免交点反向插值的舍入误差拆开同一轮廓。
    // STL 也可能用不同编号保存同坐标顶点，因此邻近九个网格桶再按实际距离合并。
    const parent=segments.map((s,i)=>i),topology=new Map(),buckets=new Map(),nodes=[],owners=[];
    function find(i){while(parent[i]!==i){parent[i]=parent[parent[i]];i=parent[i];}return i;}
    function endpoint(p,key){
      if(topology.has(key))return topology.get(key);
      const cell=p.slice(1).map((v,j)=>Math.floor((v-lo[j+1])/eps));let node;
      for(let dy=-1;dy<=1&&node===undefined;dy++)for(let dz=-1;dz<=1&&node===undefined;dz++){
        const candidates=buckets.get((cell[0]+dy)+','+(cell[1]+dz))||[];
        node=candidates.find(i=>distance(p,nodes[i])<=eps);
      }
      if(node===undefined){node=nodes.length;nodes.push(p);const bucket=cell.join(',');if(!buckets.has(bucket))buckets.set(bucket,[]);buckets.get(bucket).push(node);}
      topology.set(key,node);return node;
    }
    segments.forEach((s,i)=>{s._nodes=[s.a,s.b].map((p,j)=>endpoint(p,s._keys[j]));s._nodes.forEach(node=>{if(owners[node]!==undefined)parent[find(i)]=find(owners[node]);else owners[node]=i;});});
    const groups=new Map();segments.forEach((s,i)=>{const k=find(i);if(!groups.has(k))groups.set(k,[]);groups.get(k).push(s);});
    const components=[];
    groups.forEach(group=>{
      const perimeter=group.reduce((sum,s)=>sum+s.length,0),center=[x,0,0];
      group.forEach(s=>{center[1]+=(s.a[1]+s.b[1])*s.length/(2*perimeter);center[2]+=(s.a[2]+s.b[2])*s.length/(2*perimeter);});
      const degree=new Map();group.forEach(s=>s._nodes.forEach(node=>degree.set(node,(degree.get(node)||0)+1)));
      const closed=[...degree.values()].every(d=>d===2),branched=[...degree.values()].some(d=>d>2);
      const component=components.length,info={center,perimeter,segmentCount:group.length,closed,branched,invalidAngleSegments:0};components.push(info);
      // 模块三：周向角要求闭合无分支轮廓，且交线不能穿过角度中心或退化为径向线。
      // 空间交线、面元数值和未解状态始终保留；几何角度无效不等于物理量未求解。
      group.forEach(s=>{
        s.component=component;const ay=s.a[1]-center[1],az=s.a[2]-center[2],by=s.b[1]-center[1],bz=s.b[2]-center[2],dy=by-ay,dz=bz-az;
        const t=Math.max(0,Math.min(1,-(ay*dy+az*dz)/(s.length*s.length))),centerDistance=Math.hypot(ay+t*dy,az+t*dz),ra=Math.hypot(ay,az),rb=Math.hypot(by,bz);
        s.angleValid=false;s.angleReason=branched?'轮廓存在分支':!closed?'轮廓开放':centerDistance<=eps?'交线经过角度中心':Math.abs(ay*bz-az*by)<=eps*(ra+rb)?'交线为径向线，角度跨度退化':'';
        s.angleValid=!s.angleReason;s.thetaA=s.angleValid?angle(s.a,center):null;s.thetaB=s.angleValid?angle(s.b,center):null;if(!s.angleValid)info.invalidAngleSegments++;
        delete s._keys;delete s._nodes;
      });
    });
    return {x,segments,components,coplanarFaces,traceSide:trace===1?'+X':'-X'};
  }
  // 跨过 0°/360° 的一段拆开绘制，禁止把曲线连过整圈或跨过未解段。
  function spans(segment){
    if(segment.angleValid===false||!finite(segment.thetaA)||!finite(segment.thetaB))return [];
    let a=segment.thetaA,b=segment.thetaB,d=b-a;if(d>180)b-=360;else if(d<-180)b+=360;
    let lo=Math.min(a,b),hi=Math.max(a,b);
    if(lo<0)return [[lo+360,360],[0,hi]];
    if(hi>360)return [[lo,360],[0,hi-360]];
    return [[lo,hi]];
  }
  function sample(section,field){return section.segments.map(s=>({...s,value:finite(field.values[s.face])&&field.status[s.face]>0?field.values[s.face]:null,status:finite(field.values[s.face])?field.status[s.face]:0,method:field.methodNames[field.methodIndex[s.face]],scope:field.scopeNames?field.scopeNames[field.scopeIndex[s.face]]:field.scope}));}
  function csv(section,field,caption,condition){
    const rows=[['截面 X / m','轮廓编号','面元编号','角度A / deg','角度B / deg','A_Y / m','A_Z / m','B_Y / m','B_Z / m','交线长度 / m',field.label+' / '+field.unit,'状态','算法','适用范围','数据时刻','工况','网格环重合时取值侧','角度状态']];
    sample(section,field).forEach(s=>rows.push([section.x,s.component+1,s.face+1,s.thetaA===null?'':s.thetaA,s.thetaB===null?'':s.thetaB,s.a[1],s.a[2],s.b[1],s.b[2],s.length,s.value===null?'':s.value,s.status===1&&typeof field.computedStateLabel==='string'&&field.computedStateLabel.trim()?field.computedStateLabel:['未求解','工程计算','给定条件'][s.status],s.method,s.scope,caption,condition,section.traceSide,s.angleValid?'有效周向角':s.angleReason]));
    return '\ufeff'+rows.map(row=>row.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n');
  }
  return {geometry,spans,sample,csv};
});
