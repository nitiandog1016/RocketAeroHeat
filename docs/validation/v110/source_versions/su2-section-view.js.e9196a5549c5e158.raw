/* 模块用途：在主工作台几何画布下显示真实截面周向曲线；不计算或补全物理场。 */
(function(root,factory){const api=factory(typeof module==='object'&&module.exports?require('./ale-section-core.js'):root.RocketSection);if(typeof module==='object'&&module.exports)module.exports=api;else root.RocketSU2Section=api;})(typeof globalThis!=='undefined'?globalThis:this,function(core){
  'use strict';const finite=Number.isFinite;
  function compute(data,x){
    const f=data.faceField,view=data.sectionView;
    if(!f||!view||!finite(x))throw Error('截面需要当前逐面场和有限X。');
    if(![0,1,2].includes(view.valueStatus)||typeof view.method!=='string'||typeof view.scope!=='string')throw Error('截面状态/方法缺失。');
    const section=core.geometry(data.vertices,data.faces,x),field={label:f.name,unit:f.unit,values:f.values,
      status:f.values.map(v=>finite(v)?view.valueStatus:0),methodNames:[view.method],methodIndex:f.values.map(()=>0),
      scope:view.scope,computedStateLabel:'研究计算值；未授予物理资格'};
    return {x,revision:data.revision,fieldKey:view.fieldKey,segments:section.segments.map(s=>[s.a,s.b]),
      samples:core.sample(section,field),components:section.components,coplanarFaces:section.coplanarFaces,
      traceSide:section.traceSide,csv:core.csv(section,field,data.caption,view.condition)};
  }
  function mount(element,changed,exported){
    let data=null,state=null,face=-1,sequence=0;
    const slider=element.querySelector('[data-section-slider]'),number=element.querySelector('[data-section-x]'),
      summary=element.querySelector('[data-section-summary]'),canvas=element.querySelector('canvas');
    const fmt=v=>finite(v)?Number(v).toPrecision(5):'未知';
    function draw(){
      if(!state||element.hidden)return;
      const w=Math.max(280,canvas.clientWidth),h=Math.max(110,canvas.clientHeight),dpr=Math.min(devicePixelRatio||1,2);
      canvas.width=w*dpr;canvas.height=h*dpr;const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.clearRect(0,0,w,h);
      const l=76,r=w-28,t=23,b=h-25,range=data.faceField.range,colors=['#256686','#876c48','#567c6d'];
      ctx.font='11px system-ui';ctx.fillStyle='#526370';ctx.strokeStyle='#d1d9df';ctx.lineWidth=1;
      ctx.fillText(data.faceField.name+' / '+data.faceField.unit+' · 与云图使用同一范围',l,13);
      for(let i=0;i<=4;i++){const x=l+(r-l)*i/4;ctx.beginPath();ctx.moveTo(x,t);ctx.lineTo(x,b);ctx.stroke();ctx.fillText((i*90)+'°',x-10,b+17);}
      for(let i=0;i<=2;i++){const y=b-(b-t)*i/2;ctx.beginPath();ctx.moveTo(l,y);ctx.lineTo(r,y);ctx.stroke();ctx.fillText(fmt(range[0]+(range[1]-range[0])*i/2),2,y+4);}
      ctx.save();ctx.beginPath();ctx.rect(l,t,r-l,b-t);ctx.clip();let count=0;
      state.samples.forEach(s=>{if(!finite(s.value))return;core.spans(s).forEach(([a,z])=>{
        const y=b-(b-t)*(s.value-range[0])/(range[1]-range[0]);ctx.strokeStyle=colors[s.component%colors.length];ctx.lineWidth=s.face===face?2.8:1.4;
        ctx.beginPath();ctx.moveTo(l+(r-l)*a/360,y);ctx.lineTo(l+(r-l)*z/360,y);ctx.stroke();count++;
      });});ctx.restore();
      if(!count){ctx.fillStyle='#63717e';ctx.fillText('当前截面无可绘值/有效周向角；未知保持断开',l+10,(t+b)/2);}
      const outside=state.samples.filter(s=>finite(s.value)&&(s.value<range[0]||s.value>range[1])).length;
      summary.textContent='X '+fmt(state.x)+' m · '+state.segments.length+' 条交线 / '+state.components.length+' 个轮廓 · +Y为0°，向+Z为90° · '+(outside?outside+'段超出当前固定范围 · ':'')+'面元常值，未作曲线平滑';
    }
    function update(x){
      if(!data)return;const lo=Number(number.min),hi=Number(number.max);
      if(!finite(x)||x<lo||x>hi){number.value=state?String(state.x):'';return;}
      state=compute(data,x);state.sequence=++sequence;number.value=String(x);slider.value=String(x);
      draw();changed(state);
    }
    slider.addEventListener('input',()=>update(Number(slider.value)));
    number.addEventListener('change',()=>update(Number(number.value)));
    element.querySelector('[data-section-export]').addEventListener('click',()=>{if(state)exported(state);});
    new ResizeObserver(draw).observe(canvas);
    return {load(next){
      data=next;const enabled=Boolean(next.sectionView&&next.faceField);element.hidden=!enabled;
      if(!enabled){state=null;changed(null);return;}
      const xs=next.vertices.map(v=>v[0]),lo=Math.min(...xs),hi=Math.max(...xs);
      number.min=slider.min=String(lo);number.max=slider.max=String(hi);slider.step=String((hi-lo)/2000||1);
      face=Number.isInteger(next.sectionView.faceIndex)?next.sectionView.faceIndex:-1;
      const x=finite(next.sectionView.x)?next.sectionView.x:(state?state.x:(lo+hi)/2);
      update(Math.max(lo,Math.min(hi,x)));
    },getState:()=>state};
  }
  return {compute,mount};
});
