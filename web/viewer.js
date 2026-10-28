/* RocketAeroHeat 离线工作台：WebGL 渲染、真实三角形拾取、逐变量溯源。 */
'use strict';
// 模块零：逐变量的科学身份说明。旧JSON缺少新增元数据时明确缺失，不推断工况或借用热流原因。
(function(root,factory){const api=factory();if(typeof module==='object'&&module.exports)module.exports=api;else root.RocketViewerState=api;})(typeof globalThis!=='undefined'?globalThis:this,function(){
  const finite=x=>typeof x==='number'&&Number.isFinite(x);
  const hasHistory=field=>Boolean(field.history&&!Array.isArray(field.history)&&typeof field.history==='object');
  function state(field,i){return [1,2].includes(field.status[i])&&finite(field.values[i])?field.status[i]:0;}
  // 研究场可明确显示“未收敛数值”，旧数据继续采用原状态名称。
  function stateLabel(field,s){return s===1&&typeof field.computedStateLabel==='string'&&field.computedStateLabel.trim()?field.computedStateLabel:['未求解','工程计算','给定条件'][s]||'未求解';}
  function coverageTitle(field){return field.label+' · 面积覆盖';}
  // 面积栏遵循快照明确的研究资格；有有限数值不等于获得正式求解资格。
  function coverageLabels(model,field){
    const research=model.researchIdentity?.productionQualified===false;
    return [research?(field.computedStateLabel==='几何诊断'?'诊断数值面积':'研究数值面积'):'计算面积','给定面积','未解面积'];
  }
  function coverageNote(model,field){
    const note=stateNote(field);
    return model.researchIdentity?.productionQualified===false?note+'；本页有限研究数值不代表正式合格求解面积。':note;
  }
  function showHeatGaps(field){return field.key==='heatFlux'&&!hasHistory(field);}
  function stateNote(field){
    if(field.key==='wallSolved'&&!hasHistory(field)&&field.availability==='notProvidedInSnapshot')return '当前单点链路未执行/未提供独立壁温求解。壁温100%未解不表示对流热流全部未解。';
    if(field.stateNote)return field.stateNote;
    if(field.key==='wallInput')return '热流边界输入；给定值不是本次求解壁温。';
    if(field.key==='wallSolved'&&!hasHistory(field))return '原快照未提供独立壁温解；不能以输入壁温补齐。';
    if(field.key==='pressure'||field.key==='cp')return '表面压力的计算与给定边界分别统计；不是边界层外缘压力。';
    if(field.key==='edgePressure')return '内部外流外缘压力；局部外缘有效不代表沿程热流有效。';
    return '统计仅对应当前变量；未解原因不由其他变量推断。';
  }
  function missingReason(model,field,i){
    if(state(field,i)>0)return '';
    // 保留已记录的逐变量/实际站原因；非有限状态缺少原因时不能借用旧JSON“已计算/给定”说明。
    const invalidValue=[1,2].includes(field.status[i])&&!finite(field.values[i]);
    const invalidReason='源状态标记有效但没有有限数值；当前场量按未求解显示，不能沿用“已计算”或“给定条件”说明。';
    if(hasHistory(field)){const note=field.history.faceNote&&field.history.faceNote[i];return note||(invalidValue?invalidReason:'该实际时间站没有有效结果；不沿用单点热流拒绝原因。');}
    if(Array.isArray(field.missingReason)&&field.missingReason[i])return field.missingReason[i];
    if(invalidValue)return invalidReason;
    if(field.key==='heatFlux'&&model.missingReason&&model.missingReason[i])return model.missingReason[i];
    if(field.key==='wallSolved')return '当前快照未提供独立壁温解；不是给定壁温未知。';
    if(field.key==='wallInput')return '热流边界输入壁温未知；不由初始壁温补齐。';
    return '原结果未记录本字段具体未解原因；不以热流原因代替。';
  }
  function method(field,i){
    const name=field.methodNames[field.methodIndex[i]]||'来源未记录';
    if(state(field,i)>0||/^(未解|未求解|未计算|未接受|来源未记录)/.test(name))return name;
    return '未接受的尝试模型：'+name;
  }
  function detailTitle(field,i){
    if(state(field,i)>0)return '算法与适用范围';
    if(field.key==='wallSolved'&&!hasHistory(field)&&field.availability==='notProvidedInSnapshot')return '未执行壁温求解的说明';
    if(field.key==='wallInput')return '输入缺失与使用范围';
    return '未解原因、尝试模型与范围';
  }
  function caseCaption(model,field,index=0){
    if(hasHistory(field))return (field.history.timeNote&&field.history.timeNote[index]||'实际时间站工况未记录')+'；实际时间站，不套用单点基准工况';
    if(model.resultCondition&&model.resultCondition.caption)return model.resultCondition.caption;
    return (model.caseNote||'工况未记录')+'；旧结果缺少完整冻结工况元数据，缺失信息不推断';
  }
  // 模块用途：对照快照显式携带统一范围；未知温度没有范围，不从零值推断。
  function comparisonRange(model,key){
    const r=model.comparisonColorRanges?.[key];
    return Array.isArray(r)&&r.length===2&&r.every(finite)&&r[0]<r[1]?r.slice():null;
  }
  // 模块用途：专用研究可以明确初始变量；继续保留用户已选变量及旧快照的默认行为。
  function initialField(model,previous){
    const contains=key=>model.fields.some(f=>f.key===key);
    if(contains(previous))return previous;
    if(contains(model.defaultField))return model.defaultField;
    return contains('heatFlux')?'heatFlux':model.fields[0].key;
  }
  return {state,stateLabel,coverageTitle,coverageLabels,coverageNote,showHeatGaps,stateNote,missingReason,method,detailTitle,caseCaption,comparisonRange,initialField};
});
(function () {
  if(typeof document==='undefined')return; // Node只加载纯合同，不创建DOM或WebGL。
  const meaning=globalThis.RocketViewerState;
  const $ = id => document.getElementById(id), canvas = $('scene');
  // Viridis 为默认连续色图；保留工程蓝红色图用于与既有图对照。
  const palettes={viridis:[[68,1,84],[68,57,131],[49,104,142],[33,145,140],[53,183,121],[144,215,67],[253,231,37]],engineering:[[36,70,138],[36,127,195],[54,182,189],[147,211,155],[242,213,100],[231,119,67],[176,42,78]]};
  let colors=palettes.viridis;
  let model, field, gl, program, positions, colorBuffer, edgeBuffer, selectedBuffer;
  let center, extent, vertices, areas, centers, selected = -1, selectedPoint, bridge;
  let az = -Math.PI/2, el = .35, roll = 0, zoom = 1, pan = [0,0], bounds = [0,1], range = [0,1], drag, moved;
  let right, up, direction, sx, sy, frame = 0, baseScale = 1, pendingReport = false;
  let timeIndex=0, timeline=[], rangeHasData=false;
  let sliceX=0,sectionGeometry=null,sectionRows=[],curvePlot=null,sliceBuffer;
  const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
  const fmt=x=>Number.isFinite(x)?(Math.abs(x)>=1e5||(Math.abs(x)>0&&Math.abs(x)<.001)?x.toExponential(3):Number(x.toPrecision(5)).toString()):'未求解';
  const finite=x=>typeof x==='number'&&Number.isFinite(x);
  function problem(text){$('error').hidden=false;$('error').textContent=text;}
  // 模块一：构建 GPU 资源。没有外部库、网络资源或 MATLAB 原生三维坐标轴。
  try {
    gl=canvas.getContext('webgl',{antialias:true,preserveDrawingBuffer:true});
    if(!gl)throw Error('当前网页环境不支持 WebGL。请在支持 WebGL 的浏览器打开导出结果。');
    function shader(type,source){const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s;}
    program=gl.createProgram();
    gl.attachShader(program,shader(gl.VERTEX_SHADER,'attribute vec3 p; attribute vec3 c; uniform vec3 r,u,d; uniform vec2 scale,pan; varying vec3 color; void main(){gl_Position=vec4(dot(p,r)*scale.x+pan.x,dot(p,u)*scale.y+pan.y,-dot(p,d)*0.45,1.0);color=c;}'));
    gl.attachShader(program,shader(gl.FRAGMENT_SHADER,'precision mediump float; varying vec3 color; void main(){gl_FragColor=vec4(color,1.0);}'));
    gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(program));
    gl.useProgram(program);positions=gl.createBuffer();colorBuffer=gl.createBuffer();edgeBuffer=gl.createBuffer();selectedBuffer=gl.createBuffer();sliceBuffer=gl.createBuffer();
    gl.enable(gl.DEPTH_TEST);gl.clearColor(0,0,0,0);
  } catch(e){problem(e.message);}
  function buffer(b,values){gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(values),gl.STATIC_DRAW);}
  function attribute(name,b){const i=gl.getAttribLocation(program,name);gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.enableVertexAttribArray(i);gl.vertexAttribPointer(i,3,gl.FLOAT,false,0,0);}
  function constantColor(c){const i=gl.getAttribLocation(program,'c');gl.disableVertexAttribArray(i);gl.vertexAttrib3fv(i,c);}
  function palette(t){t=Math.max(0,Math.min(1,t))*6;const i=Math.min(5,Math.floor(t)),f=t-i;return colors[i].map((x,j)=>(x*(1-f)+colors[i+1][j]*f)/255);}
  function validate(data){
    // MATLAB 将单行矩阵/单个数编码成向量/标量，入口统一为逐面元数组。
    const list=x=>Array.isArray(x)?x:[x];
    if(data&&Array.isArray(data.faces)&&typeof data.faces[0]==='number')data.faces=[data.faces];
    if(data&&Array.isArray(data.vertices)&&typeof data.vertices[0]==='number')data.vertices=[data.vertices];
    if(data.kind!=='aerothermal-model'||data.schemaVersion!==1||!Array.isArray(data.vertices)||!data.vertices.length||!Array.isArray(data.faces)||!data.faces.length)throw Error('结果文件格式或版本无效。');
    data.fields=Array.isArray(data.fields)?data.fields:[data.fields];
    data.regions=Array.isArray(data.regions)?data.regions:[data.regions];
    data.faceArea=list(data.faceArea);data.region=list(data.region);data.missingReason=list(data.missingReason);
    data.heatingGaps=data.heatingGaps?list(data.heatingGaps):[];
    data.fields.forEach(f=>{if(f){['values','status','methodIndex','methodNames'].forEach(k=>f[k]=list(f[k]));if(f.missingReason)f.missingReason=list(f.missingReason);if(f.scopeNames&&f.scopeNames.length){f.scopeNames=list(f.scopeNames);f.scopeIndex=list(f.scopeIndex);}else{f.scopeNames=null;f.scopeIndex=null;}}});
    if(data.vertices.some(v=>v.length!==3||v.some(x=>!finite(x)))||data.faces.some(f=>f.length!==3||f.some(i=>!Number.isInteger(i)||i<0||i>=data.vertices.length)))throw Error('几何数据无效。');
    if(!data.fields.length||data.fields.some(f=>!f||f.values.length!==data.faces.length||f.status.length!==data.faces.length||f.methodIndex.length!==data.faces.length||f.status.some(s=>![0,1,2].includes(s))||f.methodIndex.some(i=>!Number.isInteger(i)||i<0||i>=f.methodNames.length)))throw Error('物理量或方法编号与面元数量不一致。');
    if(data.fields.some(f=>f.scopeNames&&(f.scopeIndex.length!==data.faces.length||f.scopeIndex.some(i=>!Number.isInteger(i)||i<0||i>=f.scopeNames.length))))throw Error('逐面元适用范围与网格不一致。');
    if(data.fields.some(f=>f.missingReason&&(f.missingReason.length!==data.faces.length||f.missingReason.some(s=>typeof s!=='string'))))throw Error('逐变量未解原因与网格不一致。');
    if(data.faceArea.length!==data.faces.length||data.faceArea.some(a=>!finite(a)||a<=0)||data.region.length!==data.faces.length||data.missingReason.length!==data.faces.length)throw Error('面积或区域信息与网格不一致。');
    let times=null;
    data.fields.forEach(f=>{
      if(!f.history||Array.isArray(f.history)&&!f.history.length){f.history=null;return;}
      const h=f.history;['time','values','status','methodIndex','methodNames','faceNote','timeNote'].forEach(k=>h[k]=list(h[k]));
      const count=data.faces.length*h.time.length;
      if(h.layout!=='time-major'||!h.time.length||h.time.some((t,j)=>!finite(t)||t<0||j>0&&t<=h.time[j-1])||h.values.length!==count||h.status.length!==count||h.methodIndex.length!==count||h.faceNote.length!==data.faces.length||h.status.some(s=>![0,1,2].includes(s))||h.methodIndex.some(i=>!Number.isInteger(i)||i<0||i>=h.methodNames.length))throw Error('热历史时间、维度或方法编号无效。');
      if(h.timeNote.length!==h.time.length)throw Error('热历史工况说明与时间轴不一致。');
      if(h.timeScope){h.timeScope=list(h.timeScope);if(h.timeScope.length!==h.time.length)throw Error('逐站适用范围与时间轴不一致。');}
      if(h.flowDirections){h.flowDirections=list(h.flowDirections);if(h.flowDirections.length!==3*h.time.length)throw Error('逐站来流方向与时间轴不一致。');}
      if(times&&JSON.stringify(times)!==JSON.stringify(h.time))throw Error('各热历史字段的时间轴不一致。');times=h.time;
      h.values=h.values.map((v,i)=>h.status[i]>0&&finite(v)?v:null);h.status=h.status.map((s,i)=>finite(h.values[i])?s:0);
    });
  }
  function load(data){
    try{validate(data);}catch(e){problem(e.message);return;}
    const previous=field&&field.key, same=model&&JSON.stringify(model.faces)===JSON.stringify(data.faces)&&JSON.stringify(model.vertices)===JSON.stringify(data.vertices);
    const oldTime=timeline[timeIndex];
    model=data;const timed=data.fields.find(f=>f.history);timeline=timed?timed.history.time:[];
    const commonOption=$('comparisonScale');if(commonOption){commonOption.disabled=!model.comparisonColorRanges;commonOption.hidden=!model.comparisonColorRanges;}
    timeIndex=same&&timeline.includes(oldTime)?timeline.indexOf(oldTime):Math.max(0,timeline.length-1);
    $('timeControls').hidden=!timeline.length;$('time').max=Math.max(0,timeline.length-1);$('time').value=timeIndex;
    let lo=[Infinity,Infinity,Infinity],hi=[-Infinity,-Infinity,-Infinity];
    data.vertices.forEach(v=>v.forEach((x,j)=>{lo[j]=Math.min(lo[j],x);hi[j]=Math.max(hi[j],x);}));
    center=lo.map((x,j)=>(x+hi[j])/2);extent=Math.max(...hi.map((x,j)=>x-lo[j]));if(!(extent>0)){problem('模型尺寸为零。');return;}
    bounds=[lo[0],hi[0]];vertices=data.vertices.map(v=>v.map((x,j)=>(x-center[j])/extent));sectionGeometry=null;
    centers=data.faces.map(f=>[0,1,2].map(j=>(data.vertices[f[0]][j]+data.vertices[f[1]][j]+data.vertices[f[2]][j])/3));areas=data.faceArea;
    if(gl){const triangles=[],edges=[];data.faces.forEach(f=>{f.forEach(i=>triangles.push(...vertices[i]));[0,1,1,2,2,0].forEach(j=>edges.push(...vertices[f[j]]));});buffer(positions,triangles);buffer(edgeBuffer,edges);}
    $('field').replaceChildren();data.fields.forEach(f=>{const o=document.createElement('option');o.value=f.key;o.textContent=f.label+' / '+f.unit;$('field').append(o);});
    $('field').value=meaning.initialField(data,previous);
    $('caseTag').textContent=data.faces.length.toLocaleString()+' 面元 · '+fmt(data.totalArea)+' m²';$('coordinates').textContent=data.coordinateNote;$('caseNote').textContent=data.caseNote;
    $('gaps').replaceChildren();data.heatingGaps.forEach(g=>{const row=document.createElement('p');row.textContent=g.reason+'：'+(100*g.areaFraction).toFixed(2)+'%（'+g.faces+' 面元）';$('gaps').append(row);});
    if(!same){selected=-1;selectedPoint=null;windwardView();zoom=1;pan=[0,0];sliceX=(bounds[0]+bounds[1])/2;}
    setSlice(sliceX);$('sliceX').min=bounds[0];$('sliceX').max=bounds[1];
    // 新对照快照默认统一色标；普通旧快照继续原选项。切换变量时重新取本变量范围。
    if(model.comparisonColorRanges)$('scale').value='comparison';
    else if($('scale').value==='comparison')$('scale').value='full';
    if(!same)fitCamera();if(gl)$('error').hidden=true;pendingReport=true;changeField();
  }
  // 模块二：线性色标、截断色标和面积覆盖率。无解值固定灰色，绝不当作零。
  // 历史字段只取实际时间站，静态字段保持原工况并明确标识。
  function atTime(f){if(!f.history)return f;const h=f.history,n=model.faces.length,start=timeIndex*n;
    return {...f,scopeNames:null,scopeIndex:null,availability:'',stateNote:'实际接受时间站的本变量；单点基准热流原因不用于解释历史缺失。',scope:f.scope+(h.timeScope?'；'+h.timeScope[timeIndex]:''),values:h.values.slice(start,start+n),status:h.status.slice(start,start+n),methodNames:h.methodNames,methodIndex:h.methodIndex.slice(start,start+n)};
  }
  function faceScope(f,i){return f.scopeNames?f.scopeNames[f.scopeIndex[i]]:f.scope;}
  function timeCaption(f){return f.history?'t = '+fmt(timeline[timeIndex])+' s':'单点工况 · 不随热响应时间变化';}
  function changeField(){if(!model)return;clearDownload('save');clearDownload('saveSurface');clearDownload('saveHistory');clearDownload('saveSection');field=atTime(model.fields.find(f=>f.key===$('field').value));$('fieldTitle').textContent=field.label+' / '+field.unit;$('fieldTime').textContent=timeCaption(field);$('scope').textContent=field.scope;
    $('caseCondition').textContent=meaning.caseCaption(model,field,timeIndex);
    $('coverageTitle').textContent=meaning.coverageTitle(field);$('coverageNote').textContent=meaning.coverageNote(model,field);
    if($('computedLabel'))$('computedLabel').textContent='● '+meaning.stateLabel(field,1);
    $('heatGapDetails').hidden=!meaning.showHeatGaps(field);
    historyScale.disabled=!field.history;if(!field.history&&$('scale').value==='history')$('scale').value='full';
    $('timeLabel').textContent=timeline.length?fmt(timeline[timeIndex])+' s · '+(timeIndex+1)+' / '+timeline.length:'';
    $('timeNote').textContent=field.history?field.history.timeNote[timeIndex]+'；实际时间站，不插值':'当前云图为单点工况，不随滑块变化';
    $('timePrev').disabled=timeIndex===0;$('timeNext').disabled=timeIndex>=timeline.length-1;
    colors=palettes[$('palette').value];document.querySelector('.gradient').style.background='linear-gradient(to top,'+colors.map(c=>'rgb('+c.join(',')+')').join(',')+')';
    const halfBand=Number($('band').value)*(bounds[1]-bounds[0])/2;
    const sectionFaces=$('scale').value==='section'&&$('sectionMode').value==='intersection'?new Set(getSection().segments.map(s=>s.face)):null;
    let values=field.values.filter((x,i)=>finite(x)&&field.status[i]>0&&($('scale').value!=='section'||(sectionFaces?sectionFaces.has(i):Math.abs(centers[i][0]-sliceX)<=halfBand))).sort((a,b)=>a-b),lo=values[0],hi=values[values.length-1];
    if($('scale').value==='history'&&field.history){lo=Infinity;hi=-Infinity;field.history.values.forEach((x,i)=>{if(finite(x)&&field.history.status[i]>0){lo=Math.min(lo,x);hi=Math.max(hi,x);}});}
    if($('scale').value==='percentile'&&values.length){lo=values[Math.floor((values.length-1)*.02)];hi=values[Math.ceil((values.length-1)*.98)];}
    if($('scale').value==='comparison'){const common=meaning.comparisonRange(model,field.key);if(common){lo=common[0];hi=common[1];}}
    if($('scale').value==='manual'){const a=$('min').value,b=$('max').value;if(a!==''&&b!==''&&Number(a)<Number(b)){lo=Number(a);hi=Number(b);}}
    rangeHasData=finite(lo)&&finite(hi);if(!rangeHasData){lo=0;hi=1;}if(lo===hi){const pad=Math.max(Math.abs(lo)*.01,.001);lo-=pad;hi+=pad;}range=[lo,hi];
    if($('scale').value!=='manual'){$('min').value=lo;$('max').value=hi;}document.querySelector('nav').classList.toggle('manual',$('scale').value==='manual');
    $('legendMin').textContent=rangeHasData?fmt(lo):'—';$('legendMax').textContent=rangeHasData?fmt(hi):'—';
    if(!rangeHasData)document.querySelector('.gradient').style.background='#4a5769';
    $('legendNote').textContent=values.length?(['percentile','section','manual'].includes($('scale').value)?'超限值截色 · 灰色未解':'灰色：无有效解'):($('scale').value==='section'?'当前截面无有效数据':field.history?'当前时间站全部未求解':'本物理量全部未求解');
    if(gl){const cs=[];field.values.forEach((v,i)=>{const c=finite(v)&&field.status[i]>0?palette((v-lo)/(hi-lo)):[.29,.34,.41];cs.push(...c,...c,...c);});buffer(colorBuffer,cs);}
    coverage();inspect();curve();historyCurve();render();
  }
  function coverage(){const f=field,totals=[0,0,0],regionMissing=Object.create(null);let total=0;
    areas.forEach((a,i)=>{const s=finite(f.values[i])?f.status[i]:0;totals[s===1?0:s===2?1:2]+=a;total+=a;if(s===0)regionMissing[model.region[i]]=(regionMissing[model.region[i]]||0)+a;});
    const pcts=totals.map(a=>a/total);$('coverageNumbers').replaceChildren();
    meaning.coverageLabels(model,f).forEach((label,i)=>{const e=document.createElement('div');e.textContent=(100*pcts[i]).toFixed(1)+'%';const s=document.createElement('small');s.textContent=label;e.append(s);$('coverageNumbers').append(e);});
    ['computed','prescribed','missing'].forEach((id,i)=>$(id).style.width=(100*pcts[i])+'%');$('regions').replaceChildren();
    model.regions.forEach(r=>{const row=document.createElement('div');row.className='region';const name=document.createElement('span');name.textContent=regionName(r.name);const metric=document.createElement('span');metric.textContent='未解 '+(100*(regionMissing[r.name]||0)/r.area).toFixed(1)+'%';row.title=fmt(r.area)+' m² · '+r.faces+' 面元';row.append(name,metric);$('regions').append(row);});
  }
  function regionName(name){return {sidewall:'侧壁',aftCap:'尾端盖',noseCap:'鼻端盖',aftTip:'尾端相邻面元',noseTip:'鼻尖相邻面元'}[name]||name;}
  // 模块三：正交相机保留真实比例；GPU 绘图与拾取使用同一相机投影。
  function camera(){direction=[Math.cos(el)*Math.cos(az),Math.cos(el)*Math.sin(az),Math.sin(el)];const br=[-Math.sin(az),Math.cos(az),0],bu=[-Math.sin(el)*Math.cos(az),-Math.sin(el)*Math.sin(az),Math.cos(el)];right=br.map((v,i)=>v*Math.cos(roll)+bu[i]*Math.sin(roll));up=bu.map((v,i)=>v*Math.cos(roll)-br[i]*Math.sin(roll));
    const aspect=canvas.clientWidth/Math.max(1,canvas.clientHeight);sx=baseScale*zoom;sy=sx*aspect;}
  function windwardView(){const current=model&&model.fields.find(f=>f.key===$('field').value),timed=current&&current.history&&model.fields.find(f=>f.history&&f.history.flowDirections),h=timed&&timed.history;const f=h?[0,1,2].map(j=>h.flowDirections[j*h.time.length+timeIndex]):model&&model.flowDirection;az=-Math.PI/2;el=.35;roll=0;if(f&&f.every(finite)&&Math.hypot(f[1],f[2])>1e-8){const y=-f[1],z=-f[2];az=y>1e-8?Math.PI/2:-Math.PI/2;el=Math.max(-1.55,Math.min(1.55,Math.atan2(z,Math.abs(y))));roll=y>1e-8?Math.PI:0;}}
  function fitCamera(){if(!vertices)return;camera();let width=0,height=0;vertices.forEach(v=>{width=Math.max(width,Math.abs(dot(v,right)));height=Math.max(height,Math.abs(dot(v,up)));});const aspect=canvas.clientWidth/Math.max(1,canvas.clientHeight);baseScale=Math.min(.79/Math.max(width,1e-8),.65/Math.max(height*aspect,1e-8));zoom=1;pan=[0,0];}
  function render(){if(frame)return;frame=requestAnimationFrame(draw);}
  function draw(){frame=0;if(!gl||!model)return;camera();const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.max(1,Math.round(canvas.clientWidth*dpr)),h=Math.max(1,Math.round(canvas.clientHeight*dpr));if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}gl.viewport(0,0,w,h);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.useProgram(program);
    [['r',right],['u',up],['d',direction]].forEach(([k,v])=>gl.uniform3fv(gl.getUniformLocation(program,k),v));gl.uniform2f(gl.getUniformLocation(program,'scale'),sx,sy);gl.uniform2fv(gl.getUniformLocation(program,'pan'),pan);
    attribute('p',positions);attribute('c',colorBuffer);gl.enable(gl.POLYGON_OFFSET_FILL);gl.polygonOffset(1,1);gl.drawArrays(gl.TRIANGLES,0,model.faces.length*3);gl.disable(gl.POLYGON_OFFSET_FILL);
    if($('wire').checked){attribute('p',edgeBuffer);constantColor([.10,.16,.23]);gl.drawArrays(gl.LINES,0,model.faces.length*6);}
    if($('showSlice').checked&&sectionGeometry){attribute('p',sliceBuffer);constantColor([.92,.94,.96]);gl.depthFunc(gl.LEQUAL);gl.drawArrays(gl.LINES,0,sectionGeometry.segments.length*2);gl.depthFunc(gl.LESS);}
    if(selected>=0){const f=model.faces[selected],lines=[];[0,1,1,2,2,0].forEach(j=>lines.push(...vertices[f[j]]));buffer(selectedBuffer,lines);attribute('p',selectedBuffer);constantColor([1,.92,.36]);gl.depthFunc(gl.LEQUAL);gl.drawArrays(gl.LINES,0,6);gl.depthFunc(gl.LESS);}
    $('axes').replaceChildren();['X','Y','Z'].forEach((t,i)=>{const s=document.createElement('span');s.textContent=t;s.style.left=(35+right[i]*30)+'px';s.style.top=(35-up[i]*30)+'px';s.style.color=['#f38c83','#84dba2','#77b9ed'][i];$('axes').append(s);});
    if(pendingReport&&bridge){pendingReport=false;const pixels=new Uint8Array(w*h*4);gl.readPixels(0,0,w,h,gl.RGBA,gl.UNSIGNED_BYTE,pixels);let visible=0;for(let i=3;i<pixels.length;i+=64)if(pixels[i])visible++;const pressure=model.fields.find(f=>f.key==='pressure');bridge.Data={kind:'rendered',faces:model.faces.length,pressureSample:pressure?pressure.values[0]:null,visiblePixelSamples:visible,webglError:gl.getError(),width:w,height:h,nonce:Date.now()};}
  }
  function project(v){return [dot(v,right)*sx+pan[0],dot(v,up)*sy+pan[1],dot(v,direction)];}
  function pick(x,y){if(!model)return;camera();const rect=canvas.getBoundingClientRect(),px=2*(x-rect.left)/rect.width-1,py=1-2*(y-rect.top)/rect.height;let nearest=-Infinity,hit=-1,weights;
    const projected=vertices.map(project);
    model.faces.forEach((f,i)=>{const a=projected[f[0]],b=projected[f[1]],c=projected[f[2]],den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]);if(Math.abs(den)<1e-14)return;
      const u=((b[1]-c[1])*(px-c[0])+(c[0]-b[0])*(py-c[1]))/den,v=((c[1]-a[1])*(px-c[0])+(a[0]-c[0])*(py-c[1]))/den,w=1-u-v;
      if(Math.min(u,v,w)<-1e-7)return;const depth=u*a[2]+v*b[2]+w*c[2];if(depth>nearest){nearest=depth;hit=i;weights=[u,v,w];}});
    if(hit<0)return;selected=hit;selectedPoint=[0,1,2].map(j=>model.faces[hit].reduce((s,i,k)=>s+model.vertices[i][j]*weights[k],0));
    setSlice(selectedPoint[0]);inspect();sliceUpdate();historyCurve();render();if(bridge)bridge.Data={kind:'selection',face:hit+1,point:selectedPoint,nonce:Date.now()};
  }
  // 模块四：显示面元常值，不伪造点插值精度；温度及其时间/来源单独列示。
  function inspect(){const box=$('probe');box.replaceChildren();if(selected<0){$('probeTitle').textContent='面元探针';$('point').textContent='点击任意面元，包括灰色区域。';return;}
    $('probeTitle').textContent='面元 #'+(selected+1)+' · '+regionName(model.region[selected]);$('point').textContent='交点 XYZ / m：'+selectedPoint.map(fmt).join(', ')+'\n面元常值；每个场量的状态和缺失原因分别列示。';
    model.fields.forEach(raw=>{const f=atTime(raw),row=document.createElement('div');row.className='probeRow';const line=document.createElement('div');line.className='valueLine';const label=document.createElement('span');label.textContent=f.label;const s=meaning.state(f,selected),value=document.createElement('span');value.className='value '+(s===2?'given':s===0?'absent':'');value.textContent=s>0?fmt(f.values[selected])+' '+f.unit:'未求解';line.append(label,value);const method=document.createElement('small');method.textContent=timeCaption(f)+' · '+meaning.stateLabel(f,s);row.append(line,method);if(s===0){const reason=document.createElement('p');reason.className='fieldReason';reason.textContent=meaning.missingReason(model,f,selected);row.append(reason);}else if(s===2){const identity=document.createElement('small');identity.textContent=meaning.stateNote(f);row.append(identity);}const detail=document.createElement('details'),summary=document.createElement('summary'),description=document.createElement('p');summary.textContent=meaning.detailTitle(f,selected);description.textContent=meaning.method(f,selected)+'；'+faceScope(f,selected)+(f.history?'；'+f.history.faceNote[selected]:'');detail.append(summary,description);row.append(detail);box.append(row);});
  }
  // 模块五：直接求 X 平面与实际三角面的交线，保留面元常值及缺口。
  function setSlice(x,keepInput=false){if(!finite(x)||!model)return;const next=Math.max(bounds[0],Math.min(bounds[1],x));if(sliceX!==next)sectionGeometry=null;sliceX=next;if(!keepInput)$('sliceX').value=Number(next.toPrecision(14));$('slice').value=bounds[1]>bounds[0]?1000*(next-bounds[0])/(bounds[1]-bounds[0]):500;}
  function getSection(){if(!sectionGeometry||sectionGeometry.x!==sliceX){sectionGeometry=RocketSection.geometry(model.vertices,model.faces,sliceX);if(gl){const lines=[];sectionGeometry.segments.forEach(s=>[s.a,s.b].forEach(v=>lines.push(...v.map((q,j)=>(q-center[j])/extent))));buffer(sliceBuffer,lines);}}return sectionGeometry;}
  function curve(){if(!model||!field)return;const geometry=getSection(),exact=$('sectionMode').value==='intersection',band=Number($('band').value)*(bounds[1]-bounds[0]);
    clearDownload('saveSection');$('curveTitle').textContent='周向分布 · '+field.label+' / '+field.unit;$('bandControl').hidden=exact;$('saveSection').disabled=!exact||!geometry.segments.length;
    sectionRows=RocketSection.sample(geometry,field);curvePlot=null;
    $('sliceLabel').textContent='X = '+fmt(sliceX)+' m'+(exact?' · '+geometry.segments.length+' 交线段':' · 带宽 '+fmt(band)+' m');
    $('curve').dataset.sectionSegments=geometry.segments.length;$('curve').dataset.sectionX=sliceX;
    const missing=sectionRows.filter(s=>s.status===0).length,angular=sectionRows.filter(s=>!s.angleValid),angleReasons=[...new Set(angular.map(s=>s.angleReason))];
    $('sectionSummary').textContent=exact?'面元常值；未解 '+missing+' 段。'+geometry.components.length+' 个轮廓，各自以周长加权中心定义角度；+Y 0°，+Z 90°。点击曲线段查看面元。'+(angular.length?' '+angular.length+' 段角度无效（'+angleReasons.join('、')+'）；空间交线与 CSV 数值保留，物理状态不变。':'')+(geometry.coplanarFaces?' 共面端盖不作为周向交线。':''):'36 个角度箱按面元面积平均，曲线缺失处断开；带内数据可能混合不同轴向位置。';
    let values=[];
    if(!exact){const origin=geometry.components.length===1?geometry.components[0].center:center,sums=Array(36).fill(0),weights=Array(36).fill(0);
      centers.forEach((c,i)=>{if(Math.abs(c[0]-sliceX)>band/2||field.status[i]===0||!finite(field.values[i]))return;const theta=(Math.atan2(c[2]-origin[2],c[1]-origin[1])+2*Math.PI)%(2*Math.PI),j=Math.min(35,Math.floor(theta/(2*Math.PI)*36));sums[j]+=field.values[i]*areas[i];weights[j]+=areas[i];});values=sums.map((s,i)=>weights[i]?s/weights[i]:null);
    }
    const valid=(exact?sectionRows.filter(s=>s.angleValid).map(s=>s.value):values).filter(finite),c=$('curve'),dpr=Math.min(window.devicePixelRatio||1,2),w=c.clientWidth,h=c.clientHeight;c.width=Math.max(1,w*dpr);c.height=Math.max(1,h*dpr);const ctx=c.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);ctx.font='10px Segoe UI';ctx.fillStyle='#8fa4bb';
    if(!valid.length){ctx.fillText(exact&&sectionRows.some(s=>finite(s.value))?'场量有效，但此截面无法定义可绘制的周向角度；查看空间交线或 CSV。':geometry.segments.length?'此截面没有该物理量的有效结果。':'此 X 平面没有非退化周向交线。',10,30);render();return;}
    let lo=Math.min(...valid),hi=Math.max(...valid);if(lo===hi){const pad=Math.max(.001,Math.abs(lo)*.01);lo-=pad;hi+=pad;}const left=64,bottom=h-19,top=9,pw=w-left-15,ph=bottom-top,px=theta=>left+theta/360*pw,py=value=>bottom-(value-lo)/(hi-lo)*ph;
    ctx.strokeStyle='#2c3a4f';ctx.lineWidth=1;[0,.5,1].forEach(t=>{const y=bottom-ph*t;ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(w-15,y);ctx.stroke();ctx.fillText(fmt(lo+(hi-lo)*t),0,y+3);});[0,90,180,270,360].forEach(v=>ctx.fillText(v+'°',px(v)-8,h-3));
    if(exact){curvePlot={left,bottom,top,pw,ph,lo,hi};sectionRows.forEach(s=>{if(s.value===null)return;ctx.strokeStyle=s.face===selected?'#e1eaf4':s.status===2?'#efbe6c':'#48d0d2';ctx.lineWidth=s.face===selected?3:1.6;ctx.setLineDash(s.component%2?[3,2]:[]);RocketSection.spans(s).forEach(([a,b])=>{ctx.beginPath();const y=py(s.value);ctx.moveTo(px(a),y);ctx.lineTo(px(b),y);ctx.stroke();});});ctx.setLineDash([]);}
    else{ctx.strokeStyle='#48d0d2';ctx.fillStyle='#48d0d2';ctx.lineWidth=1.6;let pen=false;ctx.beginPath();values.forEach((v,i)=>{if(!finite(v)){pen=false;return;}if(pen)ctx.lineTo(px((i+.5)*10),py(v));else ctx.moveTo(px((i+.5)*10),py(v));pen=true;});ctx.stroke();}
    render();
  }
  $('curve').onclick=e=>{if(!curvePlot||$('sectionMode').value!=='intersection')return;const rect=e.currentTarget.getBoundingClientRect(),x=e.clientX-rect.left,y=e.clientY-rect.top,p=curvePlot,theta=(x-p.left)/p.pw*360;let closest=null,delta=9;
    sectionRows.forEach(s=>{if(s.value===null||!RocketSection.spans(s).some(([a,b])=>theta>=a&&theta<=b))return;const sy=p.bottom-(s.value-p.lo)/(p.hi-p.lo)*p.ph,d=Math.abs(y-sy);if(d<delta){closest=s;delta=d;}});
    if(!closest)return;selected=closest.face;selectedPoint=closest.a.map((v,j)=>(v+closest.b[j])/2);inspect();curve();historyCurve();render();if(bridge)bridge.Data={kind:'selection',face:selected+1,point:selectedPoint,nonce:Date.now()};
  };
  $('saveSection').onclick=()=>{if(!model||!field||$('sectionMode').value!=='intersection')return;const text=RocketSection.csv(getSection(),field,timeCaption(field),field.history?field.history.timeNote[timeIndex]:model.caseNote);offerDownload('saveSection','RocketAeroHeat-section-X'+sliceX+'-'+field.key+'.csv','data:text/csv;charset=utf-8,'+encodeURIComponent(text),'CSV');};
  // 模块六：所选面元的离散时间曲线；有效点连线不代表重新求解或插值。
  function historyCurve(){
    const h=field&&field.history;$('historyPanel').hidden=!h;$('saveHistory').disabled=!h||selected<0;
    if(!h)return;
    const c=$('historyCurve'),w=c.clientWidth,height=c.clientHeight,dpr=Math.min(window.devicePixelRatio||1,2);
    c.width=Math.max(1,w*dpr);c.height=Math.max(1,height*dpr);const ctx=c.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,height);
    $('historyLabel').textContent=selected<0?'点击面元查看完整历史。':'面元 #'+(selected+1)+' · '+field.label+' / '+field.unit;
    if(selected<0)return;
    const n=model.faces.length,ys=h.time.map((t,j)=>h.status[j*n+selected]>0?h.values[j*n+selected]:null),valid=ys.filter(finite);
    ctx.font='10px Segoe UI';ctx.fillStyle='#91a4bb';
    if(!valid.length){ctx.fillText('该面元没有有效热历史。',8,30);return;}
    let lo=Math.min(...valid),hi=Math.max(...valid);if(lo===hi){const pad=Math.max(.001,Math.abs(lo)*.01);lo-=pad;hi+=pad;}
    const left=60,right=w-12,top=12,bottom=height-27,t0=h.time[0],span=h.time[h.time.length-1]-t0;
    const px=t=>span?left+(t-t0)/span*(right-left):(left+right)/2,py=v=>bottom-(v-lo)/(hi-lo)*(bottom-top);
    ctx.strokeStyle='#2c3a4f';[0,.5,1].forEach(f=>{const y=bottom-f*(bottom-top);ctx.beginPath();ctx.moveTo(left,y);ctx.lineTo(right,y);ctx.stroke();ctx.fillText(fmt(lo+f*(hi-lo)),0,y+3);});
    ctx.fillText(fmt(t0)+' s',left,height-7);if(span)ctx.fillText(fmt(h.time[h.time.length-1])+' s',Math.max(left,right-40),height-7);
    ctx.strokeStyle='#48d0d2';ctx.lineWidth=1.5;ctx.beginPath();let pen=false;
    ys.forEach((v,j)=>{if(!finite(v)){pen=false;return;}if(pen)ctx.lineTo(px(h.time[j]),py(v));else ctx.moveTo(px(h.time[j]),py(v));pen=true;});ctx.stroke();
    ys.forEach((v,j)=>{if(!finite(v))return;ctx.fillStyle=h.status[j*n+selected]===2?'#efbe6c':'#48d0d2';ctx.beginPath();ctx.arc(px(h.time[j]),py(v),2.5,0,Math.PI*2);ctx.fill();});
    ctx.strokeStyle='#c8d4e2';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(px(h.time[timeIndex]),top);ctx.lineTo(px(h.time[timeIndex]),bottom);ctx.stroke();
  }
  function setTime(index){if(!timeline.length)return;timeIndex=Math.max(0,Math.min(timeline.length-1,index));$('time').value=timeIndex;changeField();}
  $('time').oninput=()=>setTime(Number($('time').value));$('timePrev').onclick=()=>setTime(timeIndex-1);$('timeNext').onclick=()=>setTime(timeIndex+1);
  $('historyCurve').onclick=e=>{if(!field||!field.history)return;const r=e.currentTarget.getBoundingClientRect(),f=Math.max(0,Math.min(1,(e.clientX-r.left-60)/(r.width-72))),target=timeline[0]+f*(timeline[timeline.length-1]-timeline[0]);let closest=0;timeline.forEach((t,j)=>{if(Math.abs(t-target)<Math.abs(timeline[closest]-target))closest=j;});setTime(closest);};
  $('saveHistory').onclick=()=>{
    if(selected<0||!field.history)return;const h=field.history,n=model.faces.length,quote=v=>'"'+String(v).replace(/"/g,'""')+'"';
    const rows=[['面元编号','时间 / s',field.label+' / '+field.unit,'状态','算法','适用范围','面元说明','工况']];
    const baseScope=model.fields.find(f=>f.key===field.key).scope;
    h.time.forEach((t,j)=>{const i=j*n+selected;rows.push([selected+1,t,h.status[i]>0&&finite(h.values[i])?h.values[i]:'',meaning.stateLabel(f,h.status[i]),h.methodNames[h.methodIndex[i]],baseScope+(h.timeScope?'；'+h.timeScope[j]:''),h.faceNote[selected],h.timeNote[j]]);});
    const blob=new Blob(['\ufeff'+rows.map(r=>r.map(quote).join(',')).join('\r\n')],{type:'text/csv;charset=utf-8'});
    offerDownload('saveHistory','RocketAeroHeat-face-'+(selected+1)+'-'+field.key+'.csv',URL.createObjectURL(blob),'CSV');
  };
  // 导出保留真实下载链接，兼容不处理脚本自动下载的内嵌浏览器；不提前撤销文件。
  function offerDownload(buttonId,name,url,type,startDownload=true){let link=$(buttonId+'Download');if(!link){link=document.createElement('a');link.id=buttonId+'Download';link.className='downloadLink';$(buttonId).after(link);}if(link.href.startsWith('blob:'))URL.revokeObjectURL(link.href);link.href=url;link.download=name;link.textContent='下载 '+type;link.title=name;if(startDownload)link.click();}
  function clearDownload(buttonId){const link=$(buttonId+'Download');if(link){if(link.href.startsWith('blob:'))URL.revokeObjectURL(link.href);link.remove();}}
  const historyScale=document.createElement('option');historyScale.value='history';historyScale.textContent='全时域统一色标';$('scale').append(historyScale);
  canvas.addEventListener('pointerdown',e=>{drag={x:e.clientX,y:e.clientY,pan:e.shiftKey||e.button===2||$('panMode').checked};moved=0;canvas.setPointerCapture(e.pointerId);});
  canvas.addEventListener('pointermove',e=>{if(!drag)return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;moved+=Math.abs(dx)+Math.abs(dy);if(drag.pan){pan[0]+=2*dx/canvas.clientWidth;pan[1]-=2*dy/canvas.clientHeight;}else{az-=dx*.008;el=Math.max(-1.55,Math.min(1.55,el+dy*.008));}drag.x=e.clientX;drag.y=e.clientY;render();});
  canvas.addEventListener('pointerup',e=>{if(drag&&moved<5&&e.button===0)pick(e.clientX,e.clientY);drag=null;});canvas.addEventListener('pointercancel',()=>drag=null);canvas.addEventListener('contextmenu',e=>e.preventDefault());
  canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.1,Math.min(100,zoom*Math.exp(-e.deltaY*.001)));render();},{passive:false});
  function sliceUpdate(){clearDownload('save');if($('scale').value==='section')changeField();else curve();}
  $('field').onchange=changeField;$('scale').onchange=changeField;$('palette').onchange=changeField;$('min').onchange=changeField;$('max').onchange=changeField;$('wire').onchange=render;$('slice').oninput=()=>{setSlice(bounds[0]+Number($('slice').value)/1000*(bounds[1]-bounds[0]));sliceUpdate();};$('sliceX').onchange=()=>{if($('sliceX').value===''){setSlice(sliceX);return;}setSlice(Number($('sliceX').value));sliceUpdate();};$('band').onchange=sliceUpdate;$('sectionMode').onchange=sliceUpdate;$('showSlice').onchange=()=>{clearDownload('save');render();};
  $('sliceX').oninput=()=>{const text=$('sliceX').value,x=Number(text);if(text!==''&&finite(x)&&x>=bounds[0]&&x<=bounds[1]){setSlice(x,true);sliceUpdate();}};
  $('reset').onclick=()=>{fitCamera();render();};document.querySelectorAll('[data-view]').forEach(b=>b.onclick=()=>{roll=0;if(b.dataset.view==='windward')windwardView();else{az=b.dataset.view==='nose'?0:-Math.PI/2;el=b.dataset.view==='iso'?.6:0;if(b.dataset.view==='iso')az=-1;}fitCamera();render();});
  // 模块七：当前显示场量的完整数值交付。图内保留定量信息，配套 CSV 可核对每个面元的算法与适用边界。
  function exportOptions(){return {timeCaption:timeCaption(field),timeSeconds:field.history?timeline[timeIndex]:null,
    caseNote:meaning.caseCaption(model,field,timeIndex),resultCondition:field.history?null:model.resultCondition,range:range.slice(),rangeHasData,
    scaleLabel:$('scale').selectedOptions[0].textContent,paletteLabel:$('palette').selectedOptions[0].textContent,
    exportedAt:new Date().toISOString(),selected,selectedPoint,sceneWidth:canvas.clientWidth,sceneHeight:canvas.clientHeight,colors};}
  function surfaceDownload(options,startDownload){
    const text=RocketExport.csv(model,field,options),blob=new Blob([text],{type:'text/csv;charset=utf-8'});
    offerDownload('saveSurface',RocketExport.filename(field,options)+'.csv',URL.createObjectURL(blob),'来源 CSV',startDownload);
  }
  $('saveSurface').onclick=()=>{try{if(model&&field)surfaceDownload(exportOptions(),true);}catch(e){problem('全表面 CSV 导出失败：'+e.message);}};
  $('save').onclick=()=>{try{
    if(!gl||!model)return;draw();const options=exportOptions(),out=document.createElement('canvas'),c=out.getContext('2d');
    const plan=RocketExport.figurePlan(model,field,options,(text,size)=>{c.font=size+'px Microsoft YaHei, Segoe UI, sans-serif';return c.measureText(text).width;});
    const dpr=Math.max(1,Math.min(2,window.devicePixelRatio||1));out.width=Math.ceil(plan.width*dpr);out.height=Math.ceil(plan.height*dpr);c.scale(dpr,dpr);RocketExport.drawFigure(c,canvas,plan);
    // CSV 与 PNG 使用同一冻结快照；提供持久下载链接，避免内嵌浏览器拦截一次触发的第二个下载。
    surfaceDownload(options,false);offerDownload('save',RocketExport.filename(field,options)+'.png',out.toDataURL('image/png'),'PNG');
  }catch(e){problem('图像导出失败：'+e.message);}};
  $('openFile').onchange=e=>{const f=e.target.files[0];if(!f)return;const reader=new FileReader();reader.onload=()=>{try{load(JSON.parse(reader.result));}catch(err){problem(err.message);}};reader.readAsText(f);};
  // 模块八：标题与操作栏采用实际文档流高度；下载链接换行后色标仍位于云图区。
  const viewportResize=new ResizeObserver(()=>{
    $('viewport').style.setProperty('--viewport-header-height',Math.ceil($('viewportHeader').getBoundingClientRect().height)+'px');
    fitCamera();render();curve();historyCurve();
  });
  viewportResize.observe($('viewport'));viewportResize.observe($('viewportHeader'));
  window.setup=function(component){bridge=component;component.addEventListener('DataChanged',()=>{if(component.Data&&component.Data.kind==='aerothermal-model')load(component.Data);});if(component.Data&&component.Data.kind==='aerothermal-model')load(component.Data);component.Data={kind:'ready',nonce:Date.now()};};
  if(window.ROCKET_CASE)load(window.ROCKET_CASE);
})();
