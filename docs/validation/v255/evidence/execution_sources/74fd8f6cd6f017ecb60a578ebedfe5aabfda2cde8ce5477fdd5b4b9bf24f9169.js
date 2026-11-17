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
  // 截面 CSV 复用探针的逐变量解释；section.js 先于 viewer.js 加载，因此只在导出时解析公共模块。
  function publicMeaning(){
    return typeof module==='object'&&module.exports?require('./viewer.js'):(typeof globalThis!=='undefined'?globalThis.RocketViewerState:null);
  }
  // 研究身份只采用快照显式声明；未解行和给定条件行不冒充研究数值或正式合格结果。
  function identity(model,field,status,options){
    if(status===0)return {label:'未求解（无结果）',qualified:'不适用',note:'该行没有有效数值；不进入计算面积或正式合格面积。'};
    if(status===2)return {label:'给定条件（非求解值）',qualified:'不适用',note:field.key==='wallInput'?'输入壁温是热流边界条件，不是壁温解；不进入正式合格求解面积。':'该行是给定边界/条件，不是求解结果；不进入正式合格求解面积。'};
    const research=model&&model.researchIdentity||options&&options.researchIdentity||null;
    if(research&&research.productionQualified===false){
      const diagnostic=field.computedStateLabel==='几何诊断';
      return {label:diagnostic?'诊断数值（非正式合格）':'研究数值（非正式合格）',qualified:'否',note:diagnostic?'来源快照明确标记为几何诊断；不进入正式合格求解面积。':'来源快照明确标记为研究数值；不进入正式合格求解面积。'};
    }
    if(research&&research.productionQualified===true)return {label:'正式合格求解',qualified:'是',note:'来源快照明确标记为正式合格求解。'};
    return {label:'未声明正式资格',qualified:'未声明',note:'调用方未提供正式资格声明；不得把该数值当作正式合格结果。'};
  }
  // 兼容旧调用：只读本字段原记录；缺少原因时明确说明，不借用热流或壁温原因。
  function legacyReason(model,field,i){
    const invalid=[1,2].includes(field.status[i])&&!finite(field.values[i]);
    if(invalid)return '源状态标记有效但没有有限数值；当前场量按未求解显示，不能沿用“已计算”或“给定条件”说明。';
    if(field.history&&field.history.faceNote&&field.history.faceNote[i])return field.history.faceNote[i];
    if(Array.isArray(field.missingReason)&&field.missingReason[i])return field.missingReason[i];
    if(field.key==='heatFlux'&&model&&model.missingReason&&model.missingReason[i])return model.missingReason[i];
    if(field.key==='wallSolved')return field.availability==='notProvidedInSnapshot'?'当前单点链路未执行/未提供独立壁温解；不是给定壁温未知。':'当前快照未提供独立壁温解；不是给定壁温未知。';
    if(field.key==='wallInput')return '热流边界输入壁温未知；不由初始壁温补齐。';
    if(field.key==='edgeTemperature')return '原快照未记录外缘温度未解原因；不得借用热流原因。';
    return '原快照未记录本场量具体未解原因；不以其他变量原因代替。';
  }
  function csv(section,field,caption,condition,options={}){
    options=options||{};const model=options.model||null,meaning=publicMeaning();
    const rows=[['截面 X / m','轮廓编号','面元编号','角度A / deg','角度B / deg','A_Y / m','A_Z / m','B_Y / m','B_Z / m','交线长度 / m',field.label+' / '+field.unit,'状态','算法','适用范围','数据时刻','工况','网格环重合时取值侧','角度状态','本场量未解原因','结果身份','正式合格','身份说明']];
    sample(section,field).forEach(s=>{
      let reason='';
      if(!(s.status>0&&s.value!==null)){
        if(typeof options.missingReason==='function')reason=options.missingReason(s.face)||legacyReason(model,field,s.face);
        else if(model&&meaning&&typeof meaning.missingReason==='function')reason=meaning.missingReason(model,field,s.face);
        else reason=legacyReason(model,field,s.face);
      }
      const status=s.status===1&&typeof field.computedStateLabel==='string'&&field.computedStateLabel.trim()?field.computedStateLabel:['未求解','工程计算','给定条件'][s.status],ident=identity(model,field,s.status,options);
      rows.push([section.x,s.component+1,s.face+1,s.thetaA===null?'':s.thetaA,s.thetaB===null?'':s.thetaB,s.a[1],s.a[2],s.b[1],s.b[2],s.length,s.value===null?'':s.value,status,s.method,s.scope,caption,condition,section.traceSide,s.angleValid?'有效周向角':s.angleReason,reason,ident.label,ident.qualified,ident.note]);
    });
    return '\ufeff'+rows.map(row=>row.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n');
  }
  return {geometry,spans,sample,csv};
});
