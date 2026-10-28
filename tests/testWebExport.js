/* 单工况科研导出校核：真实模型、未解/给定区分、原精度和小画布图注。
   缺省实际箭体为19752面；第三个命令行参数可指定其他验证模型的确切面数。 */
'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),exp=require('../web/export.js');
// 独立 CSV 解码器不使用被测导出器，允许字段内逗号、引号及 CRLF。
function parseCSV(text){
  text=text.replace(/^\ufeff/,'');const rows=[];let row=[],value='',inside=false;
  for(let i=0;i<text.length;i++){
    const c=text[i];if(c==='"'){if(inside&&text[i+1]==='"'){value+='"';i++;}else inside=!inside;}
    else if(c===','&&!inside){row.push(value);value='';}
    else if((c==='\n'||c==='\r')&&!inside){if(c==='\r'&&text[i+1]==='\n')i++;row.push(value);rows.push(row);row=[];value='';}
    else value+=c;
  }
  assert.equal(inside,false,'CSV 引号必须闭合');if(value||row.length){row.push(value);rows.push(row);}return rows;
}
const vertices=[[0,0,0],[2,0,0],[0,3,0],[0,0,4]],faces=[[0,1,2],[0,1,3],[0,2,3],[1,2,3]];
const m={vertices,faces,faceArea:[3,4,6,Math.sqrt(61)],region:['nose','side','aft','肩部,"测试"\r\n第二行'],coordinateNote:'SI；鼻尖 +X；零攻角来流 -X',caseNote:'M4, α=20°；输入不是求解值',missingReason:['已计算','给定','缺失','原状态缺失']};
const f={key:'heatFlux',label:'对流热流密度',unit:'kW/m²',values:[0,1.2345678901234567,999,NaN],status:[1,2,0,1],methodIndex:[0,1,2,3],methodNames:['精确零值计算','给定,"边界"\r\n第二行','未解算法','源状态矛盾'],scope:'总范围',scopeIndex:[0,1,2,3],scopeNames:['零值范围','给定范围','未解范围','尚无有限值']};
const opts={timeSeconds:0,timeCaption:'t = 0 s',caseNote:'工况带逗号, 引号"以及\n换行',range:[0,2],rangeHasData:true,scaleLabel:'自定义色标',paletteLabel:'Viridis',exportedAt:'2026-10-13T00:00:00.000Z'};
const text=exp.csv(m,f,opts),rows=parseCSV(text),header=rows.shift(),col=name=>header.indexOf(name);
assert.ok(text.startsWith('\ufeff'));assert.equal(rows.length,4);rows.forEach(r=>assert.equal(r.length,header.length));
assert.equal(rows[0][col('数值')],'0');assert.equal(rows[0][col('状态')],'工程计算');assert.equal(rows[1][col('数值')],String(f.values[1]));assert.equal(rows[1][col('状态')],'给定条件');
assert.equal(rows[2][col('数值')],'');assert.equal(rows[2][col('状态')],'未求解');assert.equal(rows[3][col('数值')],'');assert.equal(rows[3][col('状态编号')],'0');assert.equal(rows[3][col('源状态编号')],'1');assert.match(rows[3][col('导出数据说明')],/没有有限数值/);
// 新原因列按名称读取；原40列不移位，已求零值/给定值留空，旧热流原因仍兼容。
assert.equal(col('本场量未解原因'),header.length-1);assert.equal(rows[0][col('本场量未解原因')],'');assert.equal(rows[1][col('本场量未解原因')],'');assert.equal(rows[2][col('本场量未解原因')],m.missingReason[2]);
assert.equal(rows[1][col('算法')],f.methodNames[1]);assert.equal(rows[3][col('区域')],m.region[3]);assert.equal(rows[0][col('工况')],opts.caseNote);
assert.equal(rows[0][col('时间 / s')],'0');assert.equal(rows[0][col('数据时刻类型')],'实际离散时间站');assert.equal(rows[0][col('重心 X / m')],String(2/3));assert.equal(rows[0][col('顶点3 Y / m')],'3');
assert.equal(rows[0][col('面元编号（从1开始）')],'1');assert.equal(rows[0][col('顶点1编号（从1开始）')],'1');assert.equal(rows[0][col('算法编号（从0开始）')],'0');assert.equal(rows[2][col('逐面元适用范围')],'未解范围');
const cv=exp.coverage(m,f);assert.equal(cv.computed,3);assert.equal(cv.prescribed,4);assert.equal(cv.missing,6+Math.sqrt(61));
const input={...f,key:'wallInput',label:'热流计算所用壁温',unit:'K',values:[300,310,320,330],status:[2,2,2,2]};
assert.equal(exp.coverage(m,input).computed,0);assert.equal(exp.coverage(m,input).prescribed,cv.total);
const missing={...f,key:'wallSolved',label:'求解壁温',unit:'K',values:[300,310,320,330],status:[0,0,0,0]};
assert.ok(parseCSV(exp.csv(m,missing)).slice(1).every(r=>r[col('数值')]===''&&r[col('状态')]==='未求解'));
assert.equal(parseCSV(exp.csv(m,input)).slice(1)[0][col('源热流状态说明（仅热流）')],'');
assert.throws(()=>exp.csv(m,{...f,status:[1]}),/维度/);assert.throws(()=>exp.csv(m,{...f,scopeIndex:[0,1,2,99]}),/适用范围/);
// 小画布、长中文工况、逐面元适用范围都完整换行；画布外无 fillText，原模型宽高比保持。
const measure=(s,size)=>Array.from(s).reduce((sum,c)=>sum+size*(/[\u0000-\u007f]/.test(c)?.6:1),0);
const longCase='工况包含完整定量说明。'.repeat(40),longScope='此区域属于适用范围，不得由给定值冒充计算值。'.repeat(35);
const tinyOpts={...opts,sceneWidth:270,sceneHeight:105,selected:2,selectedPoint:[0,1,2],colors:[[0,0,0],[255,255,255]],caseNote:longCase};
const plan=exp.figurePlan(m,{...f,scopeNames:['正常','给定',longScope,'缺失']},tinyOpts,measure);
assert.equal(plan.width,720);assert.ok(plan.height>tinyOpts.sceneHeight+110);assert.ok(Math.abs(plan.scene.width/plan.scene.height-270/105)<1e-12);
plan.texts.forEach(t=>{assert.ok(t.x>=0&&t.y-t.size>=0);assert.ok(t.x+measure(t.text,t.size)<=plan.width-12,JSON.stringify(t));assert.ok(t.y+4<plan.height);});
const strings=plan.texts.map(t=>t.text).join('');assert.ok(strings.includes(longCase));assert.ok(strings.includes(longScope));assert.ok(strings.includes('值：未求解；未求解'));assert.ok(!strings.includes('999 kW/m²'));
const given=exp.figurePlan(m,input,{...tinyOpts,selected:1},measure);assert.ok(given.texts.map(t=>t.text).join('').includes('值：310 K；给定条件'));
const noData=exp.figurePlan(m,missing,{...tinyOpts,rangeHasData:false,selected:1},measure);assert.ok(noData.texts.some(t=>t.text.includes('无有效数据')));assert.equal(noData.rangeHasData,false);
// 检查实际绘制使用计划尺寸、完整文本和范围；纯画布桩提供独立文字边界检查。
const painted=[],rects=[],ctx={fillRect(...r){rects.push(r);},drawImage(...r){painted.push(['image',...r]);},createLinearGradient(){return {addColorStop(){}};},strokeRect(){},fillText(s,x,y){const size=Number(this.font.match(/^\d+/)[0]);assert.ok(x+measure(s,size)<=plan.width-12);assert.ok(y+4<plan.height);painted.push(['text',s]);}};
exp.drawFigure(ctx,{},plan);assert.equal(painted.filter(x=>x[0]==='text').length,plan.texts.length);assert.deepEqual(rects[0],[0,0,plan.width,plan.height]);
// 可追溯的真实模型；缺省19,752面箭体，不使用CSV自身充当参考答案。
const actualFile=process.argv[2]||path.join(__dirname,'../outputs/v51_validation/actual_rocket/results.json');
if(fs.existsSync(actualFile)){
  const expectedFaces=process.argv[3]===undefined?19752:Number(process.argv[3]);
  assert.ok(Number.isInteger(expectedFaces)&&expectedFaces>0,'验证模型面数必须明确为正整数');
  const actual=JSON.parse(fs.readFileSync(actualFile,'utf8'));assert.equal(actual.faces.length,expectedFaces,'本次真实模型证据应保留原网格');
  for(const key of ['pressure','edgePressure','edgeTemperature','heatFlux','wallInput','wallSolved','recoveryTemperature']){
    const field={...actual.fields.find(x=>x.key===key)};
    const arows=parseCSV(exp.csv(actual,field,{timeCaption:'单点工况',caseNote:actual.caseNote})),ah=arows.shift(),ix=n=>ah.indexOf(n);
    assert.equal(arows.length,actual.faces.length);let computed=0,prescribed=0,absent=0;
    arows.forEach((r,i)=>{
      assert.equal(Number(r[ix('面元编号（从1开始）')]),i+1);assert.equal(r[ix('区域')],actual.region[i]);assert.equal(Number(r[ix('面元面积 / m²')]),actual.faceArea[i]);
      const raw=field.values[i],state=typeof raw==='number'&&Number.isFinite(raw)&&field.status[i]>0?field.status[i]:0;
      assert.equal(Number(r[ix('状态编号')]),state);assert.equal(r[ix('数值')],state?String(raw):'');assert.equal(r[ix('算法')],field.methodNames[field.methodIndex[i]]);
      if(field.missingReason)assert.equal(r[ix('本场量未解原因')],state?'':field.missingReason[i],key+'逐变量未解原因');
      assert.equal(r[ix('逐面元适用范围')],field.scopeNames&&field.scopeNames.length?field.scopeNames[field.scopeIndex[i]]:field.scope);
      if(state===1)computed+=actual.faceArea[i];else if(state===2)prescribed+=actual.faceArea[i];else absent+=actual.faceArea[i];
      const vs=actual.faces[i].map(j=>actual.vertices[j]);for(let j=0;j<3;j++){assert.equal(Number(r[ix(['重心 X / m','重心 Y / m','重心 Z / m'][j])]),(vs[0][j]+vs[1][j]+vs[2][j])/3);assert.equal(Number(r[ix('顶点'+(j+1)+'编号（从1开始）')]),actual.faces[i][j]+1);}
    });
    assert.ok(Math.abs(computed+prescribed+absent-actual.totalArea)<1e-10);if(key==='wallInput')assert.equal(computed,0);if(key==='wallSolved')assert.ok(Math.abs(absent-actual.totalArea)<1e-10);
    console.log('Actual surface CSV verified: '+key+', '+arows.length+' faces, area '+(computed+prescribed+absent)+' m².');
  }
}else if(process.argv[2])throw Error('指定的真实数据不存在：'+actualFile);
else console.log('Actual-case verification not run: pass an existing results.json as the first argument.');
console.log('Web export passed: exact full-face values/geometry/provenance, CSV escaping, zero/missing/prescribed distinction, time snapshot, adaptive small-canvas layout and full text.');
