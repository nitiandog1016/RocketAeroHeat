/* 模块一：区域结果页只消费原面元账本；不求解、不插值、不凭有限值授予资格。 */
(function(root,factory){const api=factory(typeof module==='object'&&module.exports?require('./ale-section-core.js'):root.RocketSection);if(typeof module==='object'&&module.exports)module.exports=api;else{root.RocketEvidence=api;api.mount();}})(typeof globalThis!=='undefined'?globalThis:this,function(core){
  'use strict';const finite=Number.isFinite,availabilityNames=['未知','工程计算','研究计算','给定条件','超出适用域'],qualificationNames=['未验证','验证失败','验证通过'];
  const list=v=>Array.isArray(v)?v:(v===undefined||v===null?[]:[v]);
  const at=(v,i)=>list(v)[i];
  const text=v=>v===undefined||v===null?'未记录':typeof v==='string'?v:JSON.stringify(v);
  const fmt=v=>finite(v)?Number(v).toPrecision(6):'未知';
  function fieldByKey(data,key){const f=list(data.evidenceLedger.fields).find(f=>f.key===key);if(!f)throw Error('账本缺少字段 '+key);return f;}
  function fieldRange(field){const vs=list(field.values).filter(finite);if(!vs.length)return [0,1];let lo=Infinity,hi=-Infinity;vs.forEach(v=>{lo=Math.min(lo,v);hi=Math.max(hi,v);});if(lo===hi){const gap=Math.max(Math.abs(lo)*.01,1e-6);lo-=gap;hi+=gap;}return [lo,hi];}
  function faceInfo(field,i){const a=at(field.availability,i),q=at(field.qualification,i);return {value:finite(at(field.values,i))?at(field.values,i):null,availability:a,availabilityText:availabilityNames[a]||'未知',qualification:q,qualificationText:qualificationNames[q]||'未验证',method:at(field.methodNames,at(field.methodIndex,i))||'来源未记录',scope:list(field.scopeNames).length?at(field.scopeNames,at(field.scopeIndex,i)):field.scope,reason:at(field.reason,i)||at(field.missingReason,i)||'',evidenceIDs:at(field.evidenceIDs,i)||[]};}
  function section(data,key,x){
    const f=fieldByKey(data,key),s=core.geometry(data.vertices,data.faces,x);
    const samples=s.segments.map(v=>({...v,...faceInfo(f,v.face)}));
    const rows=[['section_x_m','component_id','face_id','theta_a_deg','theta_b_deg','a_y_m','a_z_m','b_y_m','b_z_m','length_m','field_key','value','unit','availability','qualification','method','scope','reason','trace_side','angle_state']];
    samples.forEach(v=>rows.push([x,v.component+1,v.face+1,v.thetaA===null?'':v.thetaA,v.thetaB===null?'':v.thetaB,v.a[1],v.a[2],v.b[1],v.b[2],v.length,key,v.value===null?'':v.value,f.unit,v.availabilityText,v.qualificationText,v.method,v.scope,v.reason,s.traceSide,v.angleValid?'有效周向角':v.angleReason]));
    return {...s,fieldKey:key,samples,csv:'\ufeff'+rows.map(r=>r.map(v=>'"'+String(v).replace(/"/g,'""')+'"').join(',')).join('\r\n')};
  }
  function mount(){
    const $=id=>document.getElementById(id),geometrySetup=window.setup;let component,geometryListener,geometryData,data,revision=0,fieldKey='',regionName='',face=-1,position=[],ranges={},currentSection=null,geometryIdentity='',pendingFit=true;
    const geometryBridge={addEventListener:(name,fn)=>{geometryListener=fn;},get Data(){return geometryData;},set Data(message){if(!message)return;if(message.kind==='ready')return;if(message.kind==='geometry-face-probe'){inspect(message.faceIndex,message.position);return;}if(message.kind==='geometry-rendered'){send({kind:'evidence-rendered',geometry:message,...state()});return;}if(message.kind==='geometry-error')send({kind:'evidence-error',message:message.message});}};
    function send(message){if(component)component.Data={...message,revision,nonce:Date.now()+Math.random()};}
    function state(){return {fieldKey,regionName,selectedFace:face,position:[...position],ranges:{...ranges},sectionX:currentSection?currentSection.x:null,sectionInfo:currentSection};}
    function notify(){send({kind:'evidence-state',...state()});}
    function loadGeometry(){
      const f=fieldByKey(data,fieldKey);geometryData={kind:'aerothermal-geometry',schemaVersion:1,revision,vertices:data.vertices,faces:data.faces,caption:'区域结果与验证',coordinateNote:data.coordinateNote,faceField:{values:list(f.values),range:ranges[fieldKey],name:f.label,unit:f.unit,productionQualified:false},showWire:$('wire').checked,resetCamera:pendingFit,selectedFaceIDs:regionFaces(),inspectedFace:face};geometryListener();pendingFit=false;
      $('footerStatus').textContent='原面元值保留 · 数值身份与验证资格分别列示';
    }
    function regionFaces(){const r=list(data.evidenceLedger.regions).find(r=>r.name===regionName);return r?list(r.faceIDs):[];}
    function node(tag,content,cls){const n=document.createElement(tag);n.textContent=content;if(cls)n.className=cls;return n;}
    function renderTable(){
      const ledger=data.evidenceLedger,f=fieldByKey(data,fieldKey),area=f.areaSummary.totalArea;
      const box=$('coverageSummary');box.replaceChildren(node('div','全模型 '+fmt(area)+' m² · '+data.faces.length.toLocaleString()+' 个原面元','coverageLine'),node('div',f.label+'：工程 '+percent(f.areaSummary.engineeringArea,area)+' / 研究 '+percent(f.areaSummary.researchArea,area)+' / 验证通过 '+percent(f.areaSummary.qualifiedArea,area),'coverageLine'));
      const body=$('regionTable').querySelector('tbody');body.replaceChildren();
      const all={name:'全模型',area,fields:{[fieldKey]:f.areaSummary},faceIDs:[]};
      [all,...list(ledger.regions)].forEach((r,index)=>{const row=document.createElement('tr'),name=index===0?'':r.name,stats=r.fields[fieldKey];row.tabIndex=0;row.classList.toggle('selected',name===regionName);row.setAttribute('aria-selected',String(name===regionName));const first=node('td',r.name);first.append(node('small',fmt(r.area)));row.append(first);['engineeringArea','researchArea','prescribedArea','unknownArea','outOfDomainArea','qualifiedArea'].forEach(k=>{const c=node('td',percent(stats[k],r.area));c.title=fmt(stats[k])+' m²';row.append(c);});const choose=()=>selectRegion(name);row.onclick=choose;row.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose();}};body.append(row);});
      const selected=list(ledger.regions).find(r=>r.name===regionName),joint=selected?selected.joint:ledger.joint;
      $('jointSummary').textContent=(selected?selected.name:'全模型')+'：压力与热流同时有计算值 '+fmt(joint.pressureHeatFluxAvailableArea)+' m²（'+percent(joint.pressureHeatFluxAvailableArea,joint.totalArea)+'）；二者同时验证通过 '+fmt(joint.pressureHeatFluxQualifiedArea)+' m²。交集按相同原面元计算；给定值不计入。';
    }
    function percent(v,area){return area>0?(100*v/area).toFixed(2)+'%':'—';}
    function renderProbe(){
      const target=$('probeFields');target.replaceChildren();$('probeID').textContent=face<0?'未选择':'#'+(face+1);
      $('probePosition').textContent=face<0?'点击模型或选择区域定位。温度身份分别列出。':'原面元 '+(face+1)+' · '+text(at(data.region,face))+' · 面积 '+fmt(at(data.faceArea,face))+' m² · XYZ '+position.map(fmt).join(', ')+' m';
      if(face<0)return;
      // 查询默认只展开数值与两层状态；完整算法、范围和拒绝理由仍逐变量保留。
      list(data.evidenceLedger.fields).forEach(f=>{
        const p=faceInfo(f,face),block=node('div','','probeField'),head=node('div','','probeHeading');
        head.append(node('span',f.label),node('strong',fmt(p.value)+' '+f.unit));
        const status=node('div','','probeStatus');status.append(node('span',p.availabilityText,'status'),node('span',p.qualificationText,'status '+(p.qualification===1?'failed':p.qualification===2?'passed':'')));
        const details=node('details','','probeDetails');details.open=false;
        details.append(node('summary','算法与说明'),node('p','算法：'+p.method),node('p','适用范围：'+text(p.scope)));
        if(p.reason)details.append(node('p','说明：'+text(p.reason)));if(f.stateNote)details.append(node('p',f.stateNote));
        block.append(head,status,details);target.append(block);
      });
    }
    function renderEvidence(){const ledger=data.evidenceLedger,box=$('evidenceContent');box.replaceChildren(node('p','通过资格只来自与当前几何、工况、来源和原面元范围相符的证据记录；有限数值或统一方法名称不构成验证。'),node('p','接受证据 '+list(ledger.acceptedEvidence).length+' 项；拒绝证据 '+list(ledger.rejectedEvidence).length+' 项。'));const identity=node('pre',JSON.stringify(ledger.identity,null,2));box.append(identity);[['接受',ledger.acceptedEvidence],['拒绝',ledger.rejectedEvidence]].forEach(([label,items])=>list(items).forEach(item=>{const d=document.createElement('details');d.append(node('summary',label+' · '+text(item.id||item.evidenceID||item.reason||'记录')),node('pre',JSON.stringify(item,null,2)));box.append(d);}));}
    // 模块二：周向曲线直接求原面元与X平面的交线，未知断开，多轮廓保留。
    function updateSection(x){if(!data)return;const lo=Number($('sectionX').min),hi=Number($('sectionX').max);if(!finite(x)||x<lo||x>hi)return;currentSection=section(data,fieldKey,x);$('sectionX').value=$('sectionSlider').value=String(x);window.RocketGeometry.setSectionSegments(currentSection.segments.map(s=>[s.a,s.b]));drawSection();}
    function drawSection(){if(!currentSection||!data)return;const c=$('sectionCanvas'),w=Math.max(240,c.clientWidth),h=Math.max(90,c.clientHeight),dpr=Math.min(window.devicePixelRatio||1,2);c.width=w*dpr;c.height=h*dpr;const ctx=c.getContext('2d');ctx.scale(dpr,dpr);const l=76,r=w-16,t=21,b=h-23,range=ranges[fieldKey];ctx.font='10px system-ui';ctx.fillStyle='#5b6d7b';ctx.strokeStyle='#dce2e7';ctx.fillText(fieldByKey(data,fieldKey).label+' / '+fieldByKey(data,fieldKey).unit+' · 与云图同一固定范围',l,12);for(let i=0;i<=4;i++){const x=l+(r-l)*i/4;ctx.beginPath();ctx.moveTo(x,t);ctx.lineTo(x,b);ctx.stroke();ctx.fillText((90*i)+'°',x-9,b+15);}for(let i=0;i<=2;i++){const y=b-(b-t)*i/2;ctx.beginPath();ctx.moveTo(l,y);ctx.lineTo(r,y);ctx.stroke();ctx.fillText(fmt(range[0]+(range[1]-range[0])*i/2),1,y+3);}ctx.save();ctx.beginPath();ctx.rect(l,t,r-l,b-t);ctx.clip();let count=0;currentSection.samples.forEach(s=>{if(!finite(s.value))return;core.spans(s).forEach(([a,z])=>{const y=b-(b-t)*(s.value-range[0])/(range[1]-range[0]);ctx.strokeStyle=['#286883','#806c49','#527960'][s.component%3];ctx.lineWidth=s.face===face?3:1.4;ctx.beginPath();ctx.moveTo(l+(r-l)*a/360,y);ctx.lineTo(l+(r-l)*z/360,y);ctx.stroke();count++;});});ctx.restore();if(!count)ctx.fillText('无有限值或有效周向角；未知未补零',l+5,(t+b)/2);const outside=currentSection.samples.filter(s=>finite(s.value)&&(s.value<range[0]||s.value>range[1])).length;$('sectionSummary').textContent='X '+fmt(currentSection.x)+' m · '+currentSection.segments.length+' 条交线 / '+currentSection.components.length+' 轮廓 · +Y为0°，向+Z为90°'+(outside?' · '+outside+'段超出色标':'')+' · 面元常值';}
    function selectField(key){if(!data)return;const f=fieldByKey(data,key);fieldKey=key;if(!ranges[key])ranges[key]=fieldRange(f);$('fieldChoice').value=key;$('rangeLow').value=String(ranges[key][0]);$('rangeHigh').value=String(ranges[key][1]);$('fieldIdentity').textContent=f.stateNote||'数值身份与验证资格分别统计';loadGeometry();renderTable();updateSection(currentSection?currentSection.x:(Number($('sectionX').min)+Number($('sectionX').max))/2);notify();}
    function selectRegion(name){if(!data)return;if(name!==''&&!list(data.evidenceLedger.regions).some(r=>r.name===name))throw Error('未知区域');regionName=name;const ids=regionFaces();if(ids.length&&!ids.includes(face)){face=ids[0];position=data.faces[face].reduce((v,i)=>v.map((q,j)=>q+data.vertices[i][j]/3),[0,0,0]);updateSection(position[0]);}window.RocketGeometry.selectFaces(ids,face);renderTable();renderProbe();drawSection();notify();}
    function inspect(i,p){if(!data||!Number.isInteger(i)||i<0||i>=data.faces.length)throw Error('原面元编号无效');face=i;position=Array.isArray(p)&&p.length===3&&p.every(finite)?p:data.faces[i].reduce((v,k)=>v.map((q,j)=>q+data.vertices[k][j]/3),[0,0,0]);regionName=text(at(data.region,i));window.RocketGeometry.selectFaces(regionFaces(),face);renderTable();renderProbe();updateSection(position[0]);notify();}
    function load(payload){try{
      data=payload.data;revision=payload.revision;
      if(!data||!data.evidenceLedger||data.evidenceLedger.kind!=='aerothermal-evidence-ledger')throw Error('区域账本负载缺失');
      $('caseCondition').textContent=(data.resultCondition&&data.resultCondition.caption)||data.caseNote||'工况与来源见验证证据；未记录量不推断';
      const nextIdentity=JSON.stringify(data.evidenceLedger.identity.geometrySHA256);pendingFit=geometryIdentity!==nextIdentity;geometryIdentity=nextIdentity;
      const next=payload.viewState||{};ranges={...(next.ranges||{})};fieldKey=next.fieldKey||list(data.fields)[0].key;regionName=next.regionName||'';face=Number.isInteger(next.selectedFace)?next.selectedFace:-1;position=list(next.position);
      const choice=$('fieldChoice');choice.replaceChildren();list(data.fields).forEach(f=>{const o=node('option',f.label+' / '+f.unit);o.value=f.key;choice.append(o);});
      let lo=Infinity,hi=-Infinity;data.vertices.forEach(v=>{lo=Math.min(lo,v[0]);hi=Math.max(hi,v[0]);});$('sectionX').min=$('sectionSlider').min=String(lo);$('sectionX').max=$('sectionSlider').max=String(hi);$('sectionSlider').step=String((hi-lo)/2000||1);
      currentSection=null;selectField(fieldKey);if(finite(next.sectionX))updateSection(next.sectionX);renderProbe();renderEvidence();notify();
    }catch(e){$('error').hidden=false;$('error').textContent=e.message;send({kind:'evidence-error',message:e.message});}}
    // 模块三：只接受当前版本的工作台消息；导出交给MATLAB保存至用户选择的路径。
    window.setup=function(real){component=real;geometrySetup(geometryBridge);real.addEventListener('DataChanged',()=>{const p=real.Data;if(p&&p.kind==='aerothermal-evidence-view')load(p);});send({kind:'ready'});};
    $('fieldChoice').onchange=()=>selectField($('fieldChoice').value);$('applyRange').onclick=()=>{const a=Number($('rangeLow').value),b=Number($('rangeHigh').value);if(!finite(a)||!finite(b)||a>=b){$('error').hidden=false;$('error').textContent='色标须为有限且递增的范围。';return;}ranges[fieldKey]=[a,b];loadGeometry();drawSection();notify();};$('resetRange').onclick=()=>{ranges[fieldKey]=fieldRange(fieldByKey(data,fieldKey));selectField(fieldKey);};$('sectionX').onchange=()=>{updateSection(Number($('sectionX').value));notify();};$('sectionSlider').oninput=()=>{updateSection(Number($('sectionSlider').value));notify();};$('exportLedger').onclick=()=>send({kind:'evidence-export',...state()});$('exportSection').onclick=()=>send({kind:'evidence-section-export',...state()});new ResizeObserver(drawSection).observe($('sectionCanvas'));
    return {selectField,selectRegion,inspectFace:inspect,getState:state,setSectionX:x=>{updateSection(x);notify();}};
  }
  return {mount,fieldRange,faceInfo,section,availabilityNames,qualificationNames};
});
