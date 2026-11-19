/* 用工作台实际的截面、采样和 CSV 入口回读两个变量，同时核对公开版本的截面几何。 */
'use strict';
const fs = require('fs');
const path = require('path');
const vm = require('vm');
const child = require('child_process');
const assert = require('assert').strict;
const [run, output] = process.argv.slice(2).map(name => path.resolve(name));
const workbench = path.join(run, 'workbench');
const ctx = {window: {}};
vm.runInNewContext(fs.readFileSync(path.join(workbench, 'case-data.js'), 'utf8'), ctx);
const model = ctx.window.ROCKET_CASE;
const stats = JSON.parse(fs.readFileSync(path.join(run, 'audit', 'summary.json')));
const x = stats.coneRing.xM;
const section = require(path.join(workbench, 'section.js'));
const viewer = require(path.join(workbench, 'viewer.js'));
const shipped = {module: {exports: {}}};
vm.runInNewContext(child.execFileSync('git',['show','273e599d:web/section.js'],{encoding:'utf8'}),shipped);
const publishedSection = shipped.module.exports;
const cut = section.geometry(model.vertices, model.faces, x);
const reference = publishedSection.geometry(model.vertices, model.faces, x);
assert.equal(JSON.stringify(cut), JSON.stringify(reference));
assert.equal(cut.components.length, 1);
assert(cut.components[0].closed && !cut.components[0].branched && !cut.components[0].invalidAngleSegments);
assert.equal(cut.traceSide, '+X');
assert.equal(model.resultCondition.alphaDeg, 10);
assert.equal(model.researchIdentity.productionQualified, false);
assert(model.wallCoverage.fullSurfaceComputed && model.nativeEnergyBalance.belowHalfPercent);
fs.mkdirSync(output);
const fields = {};
for (const key of ['pressureGauge','heatFlux']) {
  const field = model.fields.find(f => f.key===key);
  const rows = section.sample(cut,field);
  assert.equal(JSON.stringify(rows),JSON.stringify(publishedSection.sample(reference,field)));
  const pieces = [];
  for (const r of rows) {
    assert(r.angleValid && r.status===1 && Number.isFinite(r.value));
    assert.equal(r.value,field.values[r.face]);
    assert(r.method.includes('三维半域') && viewer.stateLabel(field,r.status)==='研究数值解');
    for (const [a,b] of section.spans(r)) pieces.push({a,b,value:r.value,faceID:r.face+1});
  }
  pieces.sort((a,b)=>a.a-b.a);
  assert(Math.abs(pieces[0].a)<1e-9 && Math.abs(pieces.at(-1).b-360)<1e-9);
  for (let i=1;i<pieces.length;i++) assert(Math.abs(pieces[i].a-pieces[i-1].b)<1e-8);
  const caption=viewer.caseCaption(model,field,0);
  const csv=section.csv(cut,field,'冻结工况',caption,{model});
  assert(csv.includes('+10°') && csv.includes('研究数值解'));
  fs.writeFileSync(path.join(output,key+'.csv'),csv,{flag:'wx'});
  fields[key]={unit:field.unit,segments:rows.length,pieces,
    windward:pieces.filter(p=>p.a<=270 && p.b>=270).map(p=>p.value),
    leeward:pieces.filter(p=>p.a<=90 && p.b>=90).map(p=>p.value)};
}
fs.writeFileSync(path.join(output,'curve.json'),JSON.stringify({xM:x,traceSide:cut.traceSide,fields},null,2),{flag:'wx'});
