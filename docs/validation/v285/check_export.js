/* 用交付的离线零热流快照复查单值、手动色标和未知量的 PNG 与来源 CSV。 */
'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const dir=path.resolve(process.argv[2]),exp=require(path.join(dir,'export.js'));
const source=fs.readFileSync(path.join(dir,'case-data.js'),'utf8').trim();
const data=JSON.parse(source.slice('window.ROCKET_CASE='.length,-1));
function check(key,options,expected){
  const field=data.fields.find(f=>f.key===key);
  const plan=exp.figurePlan(data,field,options,s=>s.length*7);
  assert.deepEqual(plan.range,expected);
  const text=exp.csv(data,field,options).trimEnd().split('\r\n');
  assert.equal(text.length,5217);
  const columns=line=>Array.from(line.matchAll(/"((?:[^"]|"")*)"(?:,|$)/g),m=>m[1].replace(/""/g,'"'));
  const headers=columns(text[0]);const low=headers.indexOf('色标下限'),high=headers.indexOf('色标上限');
  assert.ok(low>=0&&high>=0);
  for(const row of text.slice(1)){
    const fields=columns(row);
    assert.equal(fields[low],options.rangeHasData?String(expected[0]):'');
    assert.equal(fields[high],options.rangeHasData?String(expected[1]):'');
  }
}
const display={rangeHasData:true,range:[-.001,.001],scaleLabel:'全模型色标',colors:[[0,0,0],[33,145,140],[255,255,255]]};
check('heatFlux',display,[0,0]);
check('heatFlux',{...display,scaleLabel:'自定义色标',range:[-2,3]},[-2,3]);
check('wallSolved',{...display,rangeHasData:false,range:[0,1]},[0,1]);
console.log('Verified 5216 rows each for uniform flux, a fixed manual scale, and unavailable solid temperature.');
