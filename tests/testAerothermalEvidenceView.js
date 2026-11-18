/* 模块用途：独立小几何核验新页面的周向原面值/资格隔离，不启动流场或MATLAB。 */
'use strict';
const assert=require('assert'),view=require('../web/evidence.js');
// 两个分离的闭合方柱轮廓，四种数据身份并存；零和负热流必须保留。
const V=[],F=[];
for(const offset of [0,4]){
  const base=V.length;
  for(const x of [0,2])for(const [y,z] of [[-1,-1],[1,-1],[1,1],[-1,1]])V.push([x,y+offset,z]);
  for(let i=0;i<4;i++){const j=(i+1)%4;F.push([base+i,base+j,base+4+j],[base+i,base+4+j,base+4+i]);}
}
const values=F.map((_,i)=>i%4===0?0:i%4===1?-2:i%4===2?null:8);
const field={key:'heatFlux',label:'热流',unit:'W/m²',values,availability:F.map((_,i)=>i%4===2?0:i%4===3?3:2),qualification:F.map((_,i)=>i%4===0?1:0),methodNames:['原算法A','给定边界'],methodIndex:F.map((_,i)=>i%4===3?1:0),scope:'局部研究',reason:F.map((_,i)=>i%4===2?'原热流未知':''),evidenceIDs:F.map(()=>[])};
const data={vertices:V,faces:F,evidenceLedger:{fields:[field]}};
const before=JSON.stringify(data),section=view.section(data,'heatFlux',1);
assert.strictEqual(section.components.length,2);assert.strictEqual(section.samples.length,F.length);
for(const sample of section.samples){
  assert.strictEqual(sample.value,values[sample.face]);assert.strictEqual(sample.availability,field.availability[sample.face]);assert.strictEqual(sample.qualification,field.qualification[sample.face]);
  assert.strictEqual(sample.a[0],1);assert.strictEqual(sample.b[0],1);assert(sample.angleValid);
}
assert(section.csv.includes('研究计算'));assert(section.csv.includes('给定条件'));assert(section.csv.includes('未知'));assert(section.csv.includes('验证失败'));assert(!section.csv.includes('验证通过'));
assert.deepStrictEqual(view.fieldRange(field),[-2,8]);assert.deepStrictEqual(view.fieldRange({values:[null,null]}),[0,1]);
assert.deepStrictEqual(view.fieldRange({values:[0,0]}),[-1e-6,1e-6]);
assert.strictEqual(JSON.stringify(data),before);assert.throws(()=>view.section(data,'absent',1));
// 模块二：用最小DOM替身核验网页状态分派；它不声称真实WebGL或鼠标验收。
class Element{
  constructor(tag=''){this.tagName=tag;this.children=[];this.clientWidth=600;this.clientHeight=140;this.classList={toggle:()=>{}};this.checked=false;}
  append(...xs){this.children.push(...xs);}replaceChildren(...xs){this.children=xs;}
  setAttribute(){}querySelector(){if(!this.body)this.body=new Element();return this.body;}
  getContext(){return new Proxy({},{get:()=>()=>{}});}
}
const elements=new Map(),element=id=>{if(!elements.has(id))elements.set(id,new Element());return elements.get(id);};
global.document={getElementById:element,createElement:tag=>new Element(tag)};global.ResizeObserver=class{observe(){}};
let geometryPayload=null,selected=null,segments=null;
global.window={devicePixelRatio:1,RocketGeometry:{getState:()=>({revision:geometryPayload?geometryPayload.revision:0}),selectFaces:(ids,face)=>{selected={ids,face};},setSectionSegments:s=>{segments=s;}}};
window.setup=bridge=>{bridge.addEventListener('DataChanged',()=>{geometryPayload=bridge.Data;selected={ids:geometryPayload.selectedFaceIDs,face:geometryPayload.inspectedFace};});bridge.Data={kind:'ready'};};
const controller=view.mount();let listener,bridgeValue;
const bridge={addEventListener:(name,fn)=>{listener=fn;},get Data(){return bridgeValue;},set Data(v){bridgeValue=v;}};
window.setup(bridge);
const stats={totalArea:16,engineeringArea:0,researchArea:8,prescribedArea:4,unknownArea:4,outOfDomainArea:0,qualifiedArea:0};
const joint={totalArea:16,pressureHeatFluxAvailableArea:8,pressureHeatFluxEngineeringArea:0,pressureHeatFluxQualifiedArea:0};
const uiData={...data,fields:[field],faceArea:F.map(()=>1),region:F.map((_,i)=>i<8?'first':'second'),evidenceLedger:{kind:'aerothermal-evidence-ledger',identity:{geometrySHA256:'same'},fields:[{...field,areaSummary:stats}],regions:['first','second'].map((name,i)=>({name,area:8,faceIDs:Array.from({length:8},(_,j)=>i*8+j),fields:{heatFlux:stats},joint})),joint,acceptedEvidence:[],rejectedEvidence:[]}};
bridgeValue={kind:'aerothermal-evidence-view',revision:1,data:uiData,viewState:{fieldKey:'heatFlux',ranges:{heatFlux:[-9,9]},sectionX:1}};listener();
assert.strictEqual(bridgeValue.kind,'evidence-state');assert.deepStrictEqual(geometryPayload.faceField.range,[-9,9]);
assert(element('jointSummary').textContent.includes('工程解 0.00000 m²（0.00%）'));
assert(element('jointSummary').textContent.includes('含研究计算的可用面积 8.00000 m²（50.00%）'));
controller.selectRegion('second');assert.strictEqual(selected.ids.length,8);assert.strictEqual(selected.face,8);
assert.deepStrictEqual(controller.getState().ranges.heatFlux,[-9,9]);assert.deepStrictEqual(geometryPayload.faceField.values,values);
controller.inspectFace(1);assert.strictEqual(controller.getState().regionName,'first');assert.strictEqual(selected.face,1);assert(segments.length>0);
assert.strictEqual(geometryPayload.faces.length,16);assert.strictEqual(geometryPayload.vertices.length,16);
const probe=element('probeFields').children[0],head=probe.children[0],status=probe.children[1],details=probe.children[2];
assert.strictEqual(head.children[0].textContent,field.label);assert.strictEqual(head.children[1].textContent,'-2.00000 W/m²');
assert.strictEqual(status.children[0].textContent,'研究计算');assert.strictEqual(status.children[1].textContent,'未验证');
assert.strictEqual(details.tagName,'details');assert.strictEqual(details.open,false);
assert(details.children.some(n=>n.textContent==='算法：原算法A'));assert(details.children.some(n=>n.textContent==='适用范围：局部研究'));
assert(!head.children.some(n=>n.tagName==='p'));assert(!status.children.some(n=>n.tagName==='p'));
const submitFace=value=>{
  element('faceNumber').value=value;
  element('faceLookup').onsubmit({preventDefault(){}});
};
for(const id of [1,F.length]){
  submitFace(String(id));
  assert.strictEqual(selected.face,id-1);assert.strictEqual(element('probeID').textContent,'#'+id);
  assert.strictEqual(controller.getState().sectionX,controller.getState().position[0]);
  assert.deepStrictEqual(controller.getState().ranges.heatFlux,[-9,9]);
}
for(const value of ['', '0','1.5','17','NaN']){
  submitFace(value);assert.strictEqual(selected.face,F.length-1);
  assert(element('lookupNote').textContent.includes('1 至 16'));
}
submitFace('1');assert.strictEqual(element('lookupNote').textContent,'');
delete joint.pressureHeatFluxEngineeringArea;
controller.selectRegion('');assert(element('jointSummary').textContent.includes('工程解 未记录'));
// 可选存档重放：通过同一个窗口的数据消息和编号表单读取实际原面，不启动浏览器。
if(process.argv[2]){
  const fs=require('node:fs'),zlib=require('node:zlib');
  const actual=JSON.parse(zlib.gunzipSync(fs.readFileSync(process.argv[2])));
  const reference=JSON.parse(fs.readFileSync(process.argv[3],'utf8'));
  bridgeValue={kind:'aerothermal-evidence-view',revision:2,data:actual,viewState:{fieldKey:'edgeTemperature'}};
  listener();
  assert.strictEqual(bridgeValue.kind,'evidence-state');
  assert.strictEqual(Number(element('faceNumber').max),reference.originalFaces);
  assert(element('jointSummary').textContent.includes('43.10%'));
  const available=['未知','工程计算','研究计算','给定条件','超出适用域'];
  for(const example of reference.probeExamples){
    submitFace(String(example.faceID));
    assert.strictEqual(selected.face,example.faceID-1);
    assert.strictEqual(element('probeID').textContent,'#'+example.faceID);
    const cards=element('probeFields').children;
    assert.strictEqual(cards.length,actual.evidenceLedger.fields.length);
    actual.evidenceLedger.fields.forEach((f,i)=>{
      const expected=example.fields[f.key],card=cards[i];
      const shown=expected.value===null?'未知':Number(expected.value).toPrecision(6);
      assert.strictEqual(card.children[0].children[1].textContent,shown+' '+expected.unit);
      assert.strictEqual(card.children[1].children[0].textContent,available[expected.availability]);
      assert.strictEqual(card.children[1].children[1].textContent,view.qualificationNames[expected.qualification]);
      assert(card.children[2].children.some(n=>n.textContent==='算法：'+expected.method));
    });
  }
  console.log(JSON.stringify({actualFaces:reference.originalFaces,queriedFaces:reference.probeExamples.map(x=>x.faceID),
    checkedFieldsPerFace:actual.evidenceLedger.fields.length,actualWebGL:false}));
}
console.log(JSON.stringify({checks:section.samples.length*6+54,components:2,originalValuesUnchanged:true,unknownNotZero:true,negativeAndZeroPreserved:true,qualificationIndependent:true,domStateDispatch:true,regionDoesNotCropOrRescale:true,compactProbeDetailsCollapsed:true,oneBasedFaceLookup:true,separateJointEngineeringArea:true,actualWebGL:false}));
