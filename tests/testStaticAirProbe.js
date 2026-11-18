/* 从实际离线页取出脚本和档案，对照原始 CSV 检查全部节点及五份内嵌导出。 */
'use strict';
const assert=require('node:assert/strict'), fs=require('node:fs'), path=require('node:path');
const crypto=require('node:crypto'), vm=require('node:vm');
const root=path.resolve(__dirname,'..'), folder=process.argv[2];
const evidence=path.join(root,'docs/validation/v261/evidence');
const page=fs.readFileSync(path.join(folder,'index.html'),'utf8');
const scripts=Array.from(page.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/g),m=>m[1]);
assert.equal(scripts.length,2);
assert.doesNotMatch(page,/<script[^>]+src\s*=|<link[^>]+href\s*=/i);
const data=JSON.parse(scripts[0]), api={exports:{}};
vm.runInNewContext(scripts[1],{module:api});
const {openArchive,probe,csv}=api.exports;
const archive=openArchive(data);
const sha=bytes=>crypto.createHash('sha256').update(bytes).digest('hex');
function csvTable(text) {
  const rows=text.replace(/^\uFEFF/,'').trimEnd().split(/\r?\n/).map(line=>{
    const cells=Array.from(line.matchAll(/(?:^|,)"((?:[^"]|"")*)"(?=,|$)/g),m=>m[1].replace(/""/g,'"'));
    assert.equal(cells.map(x=>'"'+x.replace(/"/g,'""')+'"').join(','),line);
    return cells;
  });
  const header=rows.shift();
  return rows.map(row=>{assert.equal(row.length,header.length);return Object.fromEntries(header.map((key,i)=>[key,row[i]]));});
}
let snapshots=0, wall=0, interior=0, records=0;
for (const sample of ['default','stable']) {
  for (const source of Object.values(data[sample])) {
    const bytes=fs.readFileSync(path.join(evidence,source.path));
    assert.equal(sha(bytes),source.sha256);
    assert.equal(source.text,bytes.toString('utf8'));
  }
  for (const frame of ['initial','terminal']) {
    const lines=data[sample][frame].text.trim().split(/\r?\n/);
    const headers=lines.shift().split(',');
    for (const line of lines) {
      const vals=line.split(','), raw=Object.fromEntries(headers.map((name,i)=>[name,vals[i]]));
      const id=Number(raw.id), xyz=[raw.x,raw.y,raw.z].map(Number);
      const boundary=xyz.some(x=>x===0||x===1);
      const rows=csvTable(csv(probe(archive,sample,frame,id)));
      assert.equal(rows.length,5);
      assert.deepEqual(rows.map(r=>r['物理量']),['流体节点温度','流体节点压力','给定壁温','壁面热流','求解固体壁温']);
      assert.deepEqual(rows.map(r=>r['单位']),['K','Pa','K','W/m²','K']);
      for (const row of rows) {
        assert.equal(row['配置'],sample);assert.equal(row['源帧'],frame);assert.equal(Number(row['原节点编号']),id);
        assert.deepEqual(['X_m','Y_m','Z_m'].map(k=>Number(row[k])),xyz);
        assert.equal(row['节点帧_SHA256'],data[sample][frame].sha256);
        assert.equal(row['节点帧'],sample+'/'+frame+'.csv');
        assert.equal(row['工作对象方法'],raw.method);
        assert.equal(row['工作对象记录_SHA256'],data[sample].initialization.sha256);
      }
      assert.equal(rows[0]['数值'],String(Number(raw.T_K)));assert.equal(rows[0]['方法'],raw.method);
      assert.equal(rows[1]['数值'],String(Number(raw.p_Pa)));
      assert.equal(rows[0]['变量来源'],data[sample][frame].path);
      assert.equal(rows[2]['变量来源_SHA256'],boundary ? data[sample].case.sha256 : '');
      assert.equal(rows[2]['数值'],boundary ? '300.123' : '');
      assert.equal(rows[2]['身份'],boundary ? '给定' : '不适用');
      for (const r of rows.slice(3)) {assert.equal(r['数值'],'');assert.equal(r['身份'],boundary?'未知':'不适用');}
      if (boundary) wall++; else interior++;
      snapshots++;records+=rows.length;
    }
  }
}
assert.deepEqual([snapshots,wall,interior,records],[500,392,108,2500]);
const files=JSON.parse(fs.readFileSync(path.join(folder,'report.json'),'utf8')).embeddedSelections;
for (const name of files) {
  const [,sample,frame,id]=/^static-air-(default|stable)-(initial|terminal)-node-(\d+)\.csv$/.exec(name);
  assert.equal(fs.readFileSync(path.join(folder,name),'utf8'),csv(probe(archive,sample,frame,Number(id))));
}
const old=JSON.parse(JSON.stringify(data));
old.default.initial.text=old.default.initial.text.replace(/,(?:method|TOMS748)(?=\r?$)/gm,'');
old.default.initialization.text=old.default.initialization.text.replace(/"workerMethod":"TOMS748",/,'');
const legacy=csvTable(csv(probe(openArchive(old),'default','initial',0)));
assert.equal(legacy[0]['方法'],'未记录');assert.equal(legacy[0]['工作对象方法'],'未记录');
const mismatch=JSON.parse(JSON.stringify(data));
mismatch.default.initial.text=mismatch.default.initial.text.replace('TOMS748','UNEXPECTED');
assert.throws(()=>openArchive(mismatch),/方法不一致/);
console.log(JSON.stringify({nodeSelections:snapshots,fieldRows:records,wallSelections:wall,interiorSelections:interior,
  matchedEmbeddedCSV:files.length,missingHistoricalMethod:'未记录',conflictingMethodRejected:true,externalPageDependencies:0}));
