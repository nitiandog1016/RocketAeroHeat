/* 独立几何校核：偏置方柱、面元跳变、缺失、网格环与多轮廓。 */
'use strict';
const assert=require('node:assert/strict'),section=require('../web/section.js');
function box(offset=0){
  const vertices=[[-1,1+offset,-1],[-1,3+offset,-1],[-1,3+offset,1],[-1,1+offset,1],[1,1+offset,-1],[1,3+offset,-1],[1,3+offset,1],[1,1+offset,1]];
  const faces=[[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7],[0,3,2],[0,2,1],[4,5,6],[4,6,7]];
  return {vertices,faces};
}
const m=box(),g=section.geometry(m.vertices,m.faces,.2);
assert.equal(g.segments.length,8);assert.equal(g.components.length,1);
assert.ok(Math.abs(g.components[0].perimeter-8)<1e-12);assert.ok(Math.abs(g.components[0].center[1]-2)<1e-12);
g.segments.forEach(s=>{assert.equal(s.a[0],.2);assert.equal(s.b[0],.2);assert.ok(s.length>0);});
const field={values:[10,20,30,40,50,60,70,80,90,100,110,120],status:[1,1,0,2,1,1,1,1,1,1,1,1],methodIndex:Array(12).fill(0),methodNames:['独立测试'],scope:'面元常值'};
const samples=section.sample(g,field);assert.equal(samples.find(s=>s.face===2).value,null);assert.equal(samples.find(s=>s.face===3).status,2);
assert.ok(samples.every(s=>s.value===null||s.value===field.values[s.face]));
const csv=section.csv(g,{...field,label:'压力',unit:'kPa'},'单点工况','M4');
assert.equal(csv.split('\r\n').length,9);assert.ok(csv.includes('"","未求解"'));assert.ok(csv.includes('"40","给定条件"'));assert.ok(csv.includes('"独立测试","面元常值"'));
assert.deepEqual(section.spans({thetaA:350,thetaB:10}),[[350,360],[0,10]]);
assert.deepEqual(section.spans({thetaA:10,thetaB:350}),[[350,360],[0,10]]);
assert.equal(section.geometry(m.vertices,m.faces,2).segments.length,0);
for(const x of [-1,1]){const edge=section.geometry(m.vertices,m.faces,x);assert.equal(edge.segments.length,4);assert.ok(Math.abs(edge.components[0].perimeter-8)<1e-12);assert.equal(edge.coplanarFaces,2);}
// 两个不相连的轮廓不共用错误的角度中心。
const m2=box(8),multi=section.geometry(m.vertices.concat(m2.vertices),m.faces.concat(m2.faces.map(f=>f.map(i=>i+8))),0);
assert.equal(multi.components.length,2);assert.ok(multi.components.every(c=>Math.abs(c.perimeter-8)<1e-12));
// 统一缩放和平移后仍得到相同物理截面，容差不依赖 m 单位的绝对尺寸。
const tiny=section.geometry(m.vertices.map(v=>v.map(q=>q*1e-6)),m.faces,2e-7);assert.equal(tiny.segments.length,8);
const translated=section.geometry(m.vertices.map(v=>[v[0]+100,v[1]-50,v[2]+30]),m.faces,100.2);assert.equal(translated.segments.length,8);assert.ok(Math.abs(translated.components[0].center[1]+48)<1e-10);
// 同一边正反插值的一个 ulp 误差恰好跨过旧舍入桶边界时，不得拆成两个假轮廓。
const bucketBoundary=section.geometry(m.vertices,m.faces,.000020000100000000002);
assert.equal(bucketBoundary.components.length,1);assert.equal(bucketBoundary.segments.length,8);assert.equal(bucketBoundary.components[0].closed,true);
assert.ok(Math.abs(bucketBoundary.components[0].perimeter-8)<1e-12);assert.ok(Math.abs(bucketBoundary.components[0].center[1]-2)<1e-12);assert.ok(Math.abs(bucketBoundary.components[0].center[2])<1e-12);
// 独立 STL 三角形可重复保存同坐标顶点：没有共享索引时也按真实距离连接。
const duplicateVertices=m.faces.flatMap(f=>f.map(i=>m.vertices[i])),duplicateFaces=m.faces.map((f,j)=>[3*j,3*j+1,3*j+2]);
const duplicate=section.geometry(duplicateVertices,duplicateFaces,.000020000100000000002);
assert.equal(duplicate.components.length,1);assert.equal(duplicate.components[0].closed,true);assert.ok(Math.abs(duplicate.components[0].perimeter-8)<1e-12);
function prism(polygon,xs=[-1,1]){
  const n=polygon.length,vertices=xs.flatMap(x=>polygon.map(p=>[x,...p])),faces=[];
  for(let k=0;k<xs.length-1;k++)for(let j=0;j<n;j++){const a=k*n+j,b=k*n+(j+1)%n;faces.push([a,b,b+n],[a,b+n,a+n]);}
  return {vertices,faces};
}
// 内部网格环只能采用 +X 邻侧，避免相邻两带交线重复、不同面元值混合。
const rings=prism([[-1,-1],[1,-1],[1,1],[-1,1]],[-1,0,1]),ringSection=section.geometry(rings.vertices,rings.faces,0);
assert.equal(ringSection.traceSide,'+X');assert.equal(ringSection.segments.length,4);assert.ok(ringSection.segments.every(s=>s.face>=8));assert.equal(ringSection.components[0].closed,true);
// 合法闭合 C 形截面的几何中心可位于内边上；径向交线不能冒充 180° 数值区间。
const concave=prism([[-2,-2],[2,-2],[2,-1],[0,-1],[0,1],[2,1],[2,2],[-2,2]]),concaveSection=section.geometry(concave.vertices,concave.faces,0);
assert.equal(concaveSection.components.length,1);assert.equal(concaveSection.components[0].closed,true);
assert.ok(Math.abs(concaveSection.components[0].center[1])<1e-12);assert.ok(Math.abs(concaveSection.components[0].center[2])<1e-12);
const radial=concaveSection.segments.filter(s=>!s.angleValid);assert.equal(radial.length,2);assert.ok(radial.every(s=>s.angleReason==='交线经过角度中心'&&section.spans(s).length===0));
const concaveField={values:Array(concave.faces.length).fill(7),status:Array(concave.faces.length).fill(1),methodNames:['真实测试值'],methodIndex:Array(concave.faces.length).fill(0),scope:'测试范围',label:'压力',unit:'kPa'};
assert.ok(section.sample(concaveSection,concaveField).every(s=>s.value===7&&s.status===1));
const concaveCsv=section.csv(concaveSection,concaveField,'单点','测试');assert.ok(concaveCsv.includes('"角度状态"'));assert.ok(concaveCsv.includes('"交线经过角度中心"'));
// 开放或分支轮廓保留实际交线、物理值与导出，但不提供虚构的周向角度。
const open=section.geometry([[-1,-1,0],[1,-1,0],[1,1,0]],[[0,1,2]],0);
assert.equal(open.components[0].closed,false);assert.ok(open.segments.every(s=>!s.angleValid&&s.angleReason==='轮廓开放'&&section.spans(s).length===0));
const branches=[[[0,0],[1,0]],[[0,0],[-1,0]],[[0,0],[0,1]]].map(p=>prism(p));
const branchVertices=branches.flatMap(p=>p.vertices),branchFaces=branches.flatMap((p,j)=>p.faces.slice(0,2).map(f=>f.map(i=>i+j*4)));
const branched=section.geometry(branchVertices,branchFaces,0);assert.equal(branched.components.length,1);assert.equal(branched.components[0].branched,true);assert.ok(branched.segments.every(s=>s.angleReason==='轮廓存在分支'&&section.spans(s).length===0));
console.log('Web section passed: analytical perimeter/centers, exact face values, missing/prescribed states, seam split, caps/ring traces, multiple contours, scale/translation, shared-edge rounding, duplicate vertices and radial/open/branched angle refusal.');
