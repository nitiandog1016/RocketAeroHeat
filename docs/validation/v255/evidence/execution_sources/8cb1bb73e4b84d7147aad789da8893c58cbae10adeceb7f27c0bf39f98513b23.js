/* 单工况科研导出：逐面元原值与来源 CSV，以及保持比例、自动换行的定量 PNG。 */
(function(root,factory){const api=factory(root);if(typeof module==='object'&&module.exports)module.exports=api;else root.RocketExport=api;})(typeof globalThis!=='undefined'?globalThis:this,function(root){
  'use strict';
  // 逐变量未解原因复用探针的同一解析；网页在调用时读取，兼容 export.js 先于 viewer.js 加载。
  const fieldMeaning=typeof module==='object'&&module.exports?require('./viewer.js'):null;
  const finite=v=>typeof v==='number'&&Number.isFinite(v),states=['未求解','工程计算','给定条件'];
  const format=v=>finite(v)?(Math.abs(v)>=1e5||(Math.abs(v)>0&&Math.abs(v)<.001)?v.toExponential(3):Number(v.toPrecision(5)).toString()):'未求解';
  // 模块一：同一显示时刻的数据快照。状态优先于数值，未解占位数不能进入 CSV、探针或图注。
  function faceState(field,i){return [1,2].includes(field.status[i])&&finite(field.values[i])?field.status[i]:0;}
  function stateLabel(field,s){const meaning=fieldMeaning||root.RocketViewerState;return meaning&&meaning.stateLabel?meaning.stateLabel(field,s):states[s];}
  function faceScope(field,i){return field.scopeNames&&field.scopeNames.length&&field.scopeIndex?field.scopeNames[field.scopeIndex[i]]:field.scope||'';}
  function faceMethod(field,i){return field.methodNames[field.methodIndex[i]]||'';}
  function faceMissingReason(model,field,i){
    const meaning=fieldMeaning||root.RocketViewerState;
    if(!meaning||typeof meaning.missingReason!=='function')throw Error('逐变量未解原因模块尚未加载。');
    return meaning.missingReason(model,field,i);
  }
  function check(model,field){
    if(!model||!field||!Array.isArray(model.faces)||!model.faces.length)throw Error('没有可导出的表面结果。');
    const n=model.faces.length;
    for(const key of ['values','status','methodIndex'])if(!Array.isArray(field[key])||field[key].length!==n)throw Error('导出场量维度与网格不一致：'+key);
    if(!Array.isArray(field.methodNames)||model.faceArea.length!==n||model.region.length!==n)throw Error('导出面积、区域或算法信息不完整。');
    if(field.status.some(s=>![0,1,2].includes(s))||field.methodIndex.some(i=>!Number.isInteger(i)||i<0||i>=field.methodNames.length))throw Error('导出状态或算法编号无效。');
    if(field.scopeNames&&field.scopeNames.length&&(!field.scopeIndex||field.scopeIndex.length!==n||field.scopeIndex.some(i=>!Number.isInteger(i)||i<0||i>=field.scopeNames.length)))throw Error('逐面元适用范围不完整。');
    if(field.missingReason&&(!Array.isArray(field.missingReason)||field.missingReason.length!==n||field.missingReason.some(reason=>typeof reason!=='string')))throw Error('逐变量未解原因与网格不一致。');
  }
  function coverage(model,field){
    const area=[0,0,0];model.faces.forEach((f,i)=>area[faceState(field,i)]+=model.faceArea[i]);
    const total=area.reduce((a,b)=>a+b,0);return {total,missing:area[0],computed:area[1],prescribed:area[2]};
  }
  const headers=['面元编号（从1开始）','顶点1编号（从1开始）','顶点2编号（从1开始）','顶点3编号（从1开始）',
    '顶点1 X / m','顶点1 Y / m','顶点1 Z / m','顶点2 X / m','顶点2 Y / m','顶点2 Z / m','顶点3 X / m','顶点3 Y / m','顶点3 Z / m',
    '重心 X / m','重心 Y / m','重心 Z / m','面元面积 / m²','区域','场量标识','物理量','数值','单位','状态编号','状态','源状态编号',
    '算法编号（从0开始）','算法','逐面元适用范围','源热流状态说明（仅热流）','导出数据说明','数据时刻类型','时间 / s','时刻说明','工况',
    '坐标约定','色标下限','色标上限','色标方式','色图','导出时间（UTC）','本场量未解原因'];
  function row(model,field,i,options={}){
    const f=model.faces[i],vs=f.map(j=>model.vertices[j]),center=[0,1,2].map(j=>vs.reduce((s,v)=>s+v[j],0)/3),state=faceState(field,i);
    return [i+1,...f.map(j=>j+1),...vs.flat(),...center,model.faceArea[i],model.region[i],field.key,field.label,
      state?field.values[i]:'',field.unit,state,stateLabel(field,state),field.status[i],field.methodIndex[i],faceMethod(field,i),faceScope(field,i),
      field.key==='heatFlux'&&model.missingReason?model.missingReason[i]:'',
      field.status[i]>0&&state===0?'源状态有效但没有有限数值；导出为未求解':'面元常值；不插值；给定条件不计为工程计算',
      finite(options.timeSeconds)?'实际离散时间站':'单点工况',finite(options.timeSeconds)?options.timeSeconds:'',options.timeCaption||'单点工况',
      options.caseNote===undefined?model.caseNote||'':options.caseNote,model.coordinateNote||'',
      options.rangeHasData&&options.range?options.range[0]:'',options.rangeHasData&&options.range?options.range[1]:'',
      options.scaleLabel||'',options.paletteLabel||'',options.exportedAt||'',faceMissingReason(model,field,i)];
  }
  // 模块二：RFC 4180 引号与 CRLF；保留浮点原精度、中文 BOM、缺失空字段和真实 MATLAB 面元编号。
  function csv(model,field,options={}){
    check(model,field);const quote=v=>'"'+String(v===undefined||v===null?'':v).replace(/"/g,'""')+'"',lines=[headers.map(quote).join(',')];
    model.faces.forEach((f,i)=>lines.push(row(model,field,i,options).map(quote).join(',')));
    return '\ufeff'+lines.join('\r\n')+'\r\n';
  }
  function filename(field,options={}){return 'RocketAeroHeat-'+String(field.key).replace(/[^a-zA-Z0-9_-]/g,'_')+(finite(options.timeSeconds)?'-t'+options.timeSeconds+'s':'-steady');}
  // 模块三：按实际字体宽度逐字换行；不剪裁工况、单位、来源或所选面元状态。
  function wrap(text,maxWidth,measure){
    const lines=[];
    String(text||'').split(/\r\n|\n|\r/).forEach(paragraph=>{
      let line='';for(const c of Array.from(paragraph)){if(line&&measure(line+c)>maxWidth){lines.push(line);line=c;}else line+=c;}lines.push(line);
    });return lines;
  }
  function figurePlan(model,field,options,measure){
    check(model,field);
    const width=Math.max(720,Math.round(options.sceneWidth||720)),margin=24,fontSize=12,lineHeight=19,texts=[];
    const fit=(text,x,y,maxWidth,size=fontSize,color='#b7c9d9')=>{
      const lines=wrap(text,maxWidth,t=>measure(t,size));lines.forEach((text,i)=>texts.push({text,x,y:y+i*(size===16?24:lineHeight),size,color}));
      return y+lines.length*(size===16?24:lineHeight);
    };
    let y=fit(field.label+' / '+field.unit+' · '+(options.timeCaption||'单点工况'),margin,29,width-2*margin,16,'#e1eaf4');
    y=fit(model.coordinateNote||'SI · 真实几何比例',margin,y+3,width-2*margin);
    const top=y+10,sceneWidth=width-2*margin-154,sceneHeight=Math.max(240,Math.round(options.sceneHeight||360));
    const sourceW=Math.max(1,options.sceneWidth||sceneWidth),sourceH=Math.max(1,options.sceneHeight||sceneHeight),factor=Math.min(sceneWidth/sourceW,sceneHeight/sourceH);
    const scene={x:margin+(sceneWidth-sourceW*factor)/2,y:top+(sceneHeight-sourceH*factor)/2,width:sourceW*factor,height:sourceH*factor};
    const legend={x:width-margin-111,y:top+30,width:17,height:Math.max(130,sceneHeight-90)};
    fit(field.unit,legend.x,top+8,111,12,'#e1eaf4');
    const ticks=[1,.75,.5,.25,0].map(t=>({text:options.rangeHasData?format(options.range[0]+t*(options.range[1]-options.range[0])):'—',x:legend.x+25,y:legend.y+(1-t)*legend.height+4,size:11,color:'#e1eaf4'}));texts.push(...ticks);
    fit('灰色：未求解',legend.x,legend.y+legend.height+25,111,11);
    y=top+sceneHeight+20;
    y=fit('色标：'+(options.scaleLabel||'全模型色标')+'；范围 '+(options.rangeHasData?format(options.range[0])+' ～ '+format(options.range[1])+' '+field.unit:'无有效数据')+'；'+(options.paletteLabel||'Viridis'),margin,y,width-2*margin);
    if(options.scaleLabel&&/分位/.test(options.scaleLabel))y=fit('超出分位范围的有效数值显示端点颜色；原值保留在 CSV。',margin,y,width-2*margin);
    const c=coverage(model,field),pct=a=>(100*a/c.total).toFixed(2)+'%';
    y=fit('按面积统计：'+stateLabel(field,1)+' '+pct(c.computed)+'；给定条件 '+pct(c.prescribed)+'；未求解 '+pct(c.missing)+'。',margin,y,width-2*margin);
    y=fit('工况：'+(options.caseNote===undefined?model.caseNote||'未提供':options.caseNote),margin,y,width-2*margin);
    const methods=[...new Set(model.faces.map((f,i)=>faceMethod(field,i)))];
    y=fit('算法来源共 '+methods.length+' 类；全表面逐面元数值、状态、算法和适用范围见 '+filename(field,options)+'.csv。',margin,y,width-2*margin);
    if(Number.isInteger(options.selected)&&options.selected>=0&&options.selected<model.faces.length){
      const i=options.selected,state=faceState(field,i),point=options.selectedPoint||[0,1,2].map(j=>model.faces[i].reduce((s,k)=>s+model.vertices[k][j],0)/3);
      y=fit('面元 #'+(i+1)+'；XYZ / m：'+point.map(format).join(', ')+'；值：'+(state?format(field.values[i])+' '+field.unit:'未求解')+'；'+stateLabel(field,state),margin,y+5,width-2*margin,12,'#e1eaf4');
      y=fit('该面元算法：'+faceMethod(field,i),margin,y,width-2*margin);
      y=fit('该面元适用范围：'+faceScope(field,i),margin,y,width-2*margin);
    }
    return {width,height:Math.ceil(y+margin),texts,scene,legend,rangeHasData:options.rangeHasData,colors:options.colors};
  }
  // 模块四：图像内容与颜色使用当前显示场量；图中模型保持原比例，图注和色标放在独立留白内。
  function drawFigure(ctx,source,plan){
    ctx.fillStyle='#162234';ctx.fillRect(0,0,plan.width,plan.height);const s=plan.scene;ctx.drawImage(source,s.x,s.y,s.width,s.height);
    const l=plan.legend;
    if(plan.rangeHasData){const g=ctx.createLinearGradient(0,l.y+l.height,0,l.y);plan.colors.forEach((rgb,i)=>g.addColorStop(i/(plan.colors.length-1),'rgb('+rgb.join(',')+')'));ctx.fillStyle=g;}else ctx.fillStyle='#4a5769';
    ctx.fillRect(l.x,l.y,l.width,l.height);ctx.strokeStyle='#8192a6';ctx.lineWidth=.5;ctx.strokeRect(l.x,l.y,l.width,l.height);
    ctx.textBaseline='alphabetic';plan.texts.forEach(t=>{ctx.font=t.size+'px Microsoft YaHei, Segoe UI, sans-serif';ctx.fillStyle=t.color;ctx.fillText(t.text,t.x,t.y);});
  }
  return {headers,faceState,faceScope,faceMethod,faceMissingReason,coverage,row,csv,filename,wrap,figurePlan,drawFigure};
});
