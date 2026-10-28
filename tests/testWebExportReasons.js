/* 当前单工况 CSV 的逐变量未解原因：独立解码、原列兼容及探针/导出同源。
   不求解物理场；真实结果与 v57 六面来源 fixture 均只读核对。 */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const exp=require('../web/export.js'),meaning=require('../web/viewer.js');
function parseCSV(text){
  const rows=[];let row=[],value='',quoted=false;
  text=text.replace(/^\ufeff/,'');
  for(let i=0;i<text.length;i++){
    const c=text[i];
    if(c==='"'){if(quoted&&text[i+1]==='"'){value+='"';i++;}else quoted=!quoted;}
    else if(c===','&&!quoted){row.push(value);value='';}
    else if((c==='\r'||c==='\n')&&!quoted){if(c==='\r'&&text[i+1]==='\n')i++;row.push(value);rows.push(row);row=[];value='';}
    else value+=c;
  }
  assert.equal(quoted,false);if(value||row.length){row.push(value);rows.push(row);}return rows;
}
const m={vertices:[[0,0,0],[1,0,0],[0,1,0]],faces:[[0,1,2],[0,1,2],[0,1,2]],faceArea:[.5,.5,.5],region:['鼻部','侧壁','底面'],missingReason:['旧热流A','旧热流B','旧热流C'],caseNote:'合成单工况',coordinateNote:'SI'};
const f={key:'heatFlux',label:'对流热流密度',unit:'W/m²',values:[0,0,null],status:[1,2,0],methodIndex:[0,0,0],methodNames:['来源'],scope:'范围',missingReason:['不得导出已算原因','不得导出给定原因','本面未知，不能借用上游\r\n第二行 "原因"']};
const oldHeaders=exp.headers.slice(0,-1);
assert.equal(oldHeaders.length,40,'原有40列及顺序必须保留');assert.equal(exp.headers.at(-1),'本场量未解原因');
let decoded=parseCSV(exp.csv(m,f)),header=decoded.shift(),reason=header.indexOf('本场量未解原因');
assert.deepEqual(header.slice(0,-1),oldHeaders);assert.equal(reason,40);
assert.equal(decoded[0][reason],'');assert.equal(decoded[1][reason],'');assert.equal(decoded[2][reason],f.missingReason[2]);
assert.equal(decoded[2][header.indexOf('源热流状态说明（仅热流）')],m.missingReason[2],'旧原始热流列不改写');
assert.equal(exp.faceMissingReason(m,{...f,missingReason:undefined},2),'旧热流C');
const edge={...f,key:'edgeTemperature',missingReason:undefined};
assert.match(exp.faceMissingReason(m,edge,2),/本字段具体未解原因/);assert.ok(!exp.faceMissingReason(m,edge,2).includes('旧热流C'));
assert.match(exp.faceMissingReason(m,{...edge,key:'wallInput'},2),/输入壁温未知/);
assert.match(exp.faceMissingReason(m,{...edge,key:'wallSolved'},2),/独立壁温解/);
assert.equal(exp.faceMissingReason(m,{...f,history:[]},2),f.missingReason[2],'空history数组仍是单点字段');
const timed={...f,history:{faceNote:['实际站A','实际站B','实际站原因，不借静态原因']}};
assert.equal(exp.faceMissingReason(m,timed,2),'实际站原因，不借静态原因');
assert.match(exp.faceMissingReason(m,{...timed,history:{faceNote:[]}},2),/实际时间站/);
assert.throws(()=>exp.csv(m,{...f,missingReason:['短数组']}),/逐变量未解原因/);
assert.throws(()=>exp.csv(m,{...f,missingReason:['','',{}]}),/逐变量未解原因/);
const contradiction={...edge,values:[0,0,NaN],status:[1,2,1],missingReason:['','','非有限诊断值，没有接受解']};
decoded=parseCSV(exp.csv(m,contradiction));header=decoded.shift();assert.equal(decoded[2][header.indexOf('数值')],'');
assert.equal(decoded[2][reason],'非有限诊断值，没有接受解');assert.match(decoded[2][header.indexOf('导出数据说明')],/没有有限数值/);
const legacyInvalid={...contradiction,key:'heatFlux',missingReason:undefined};
const legacyRow=parseCSV(exp.csv({...m,missingReason:['已计算','给定条件','已计算']},legacyInvalid))[3];
assert.match(legacyRow[reason],/按未求解显示/);assert.notEqual(legacyRow[reason],'已计算');
// 浏览器脚本加载顺序的只读合同：export先加载、viewer后加载时也通过同一个meaning函数。
const vm=require('node:vm'),browser={};browser.globalThis=browser;
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../web/export.js'),'utf8'),browser);
vm.runInNewContext(fs.readFileSync(path.join(__dirname,'../web/viewer.js'),'utf8'),browser);
assert.equal(browser.RocketExport.faceMissingReason(m,f,2),browser.RocketViewerState.missingReason(m,f,2));
assert.ok(fs.readFileSync(path.join(__dirname,'../web/viewer.js'),'utf8').includes('reason.textContent=meaning.missingReason(model,f,selected)'),'探针继续使用同一个纯解析函数');
const files=process.argv.slice(2),summaries=[];
for(const file of files){
  const bytes=fs.readFileSync(file),model=JSON.parse(bytes),keys=['pressure','edgePressure','edgeTemperature','heatFlux','wallInput','wallSolved','recoveryTemperature'];
  const fields=[];
  for(const key of keys){
    const field=model.fields.find(x=>x.key===key);assert.ok(field,key);
    const rows=parseCSV(exp.csv(model,field));const h=rows.shift(),ri=h.indexOf('本场量未解原因'),si=h.indexOf('状态编号'),vi=h.indexOf('数值');
    assert.equal(rows.length,model.faces.length);assert.equal(ri,h.length-1);
    let missing=0;
    rows.forEach((row,i)=>{
      assert.equal(row.length,h.length);
      const valid=[1,2].includes(field.status[i])&&typeof field.values[i]==='number'&&Number.isFinite(field.values[i]);
      const expected=valid?'':field.missingReason[i];
      assert.equal(row[ri],expected,key+' face '+(i+1));
      assert.equal(row[ri],meaning.missingReason(model,field,i),key+' probe/CSV '+(i+1));
      if(valid){assert.equal(row[si],String(field.status[i]));assert.equal(row[vi],String(field.values[i]));}
      else{missing++;assert.equal(row[si],'0');assert.equal(row[vi],'');assert.ok(row[ri].length>0);}
    });
    fields.push({key,rows:rows.length,missingReasonRows:missing});
    console.log('CSV_REASON_FIELD_PASS '+path.basename(path.dirname(file))+' '+key+' rows='+rows.length+' missing='+missing);
  }
  summaries.push({file,sha256:crypto.createHash('sha256').update(bytes).digest('hex'),faces:model.faces.length,fields});
}
if(process.env.ROCKET_UI_REASON_REPORT)fs.writeFileSync(process.env.ROCKET_UI_REASON_REPORT,JSON.stringify({status:'passed',syntheticCases:['computed-zero','prescribed-zero','per-field-reason','legacy-heat-only-fallback','independent-temperature-fallback','empty-history','actual-station-precedence','contradictory-nonfinite','invalid-reason-schema','browser-load-order'],sourceFixtures:summaries,physicalRecomputation:false,newSolvedSurfaceAreaM2:0},null,2)+'\n');
console.log('WEB_EXPORT_VARIABLE_REASON_PASS: original 40 columns preserved; appended current-field reason matches probe; valid zero/prescribed blank; legacy/history/schema and actual snapshots verified.');
