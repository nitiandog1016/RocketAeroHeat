/* 模块用途：同窗只读研究页桥接，按当前模型显示几何/网格/已保存条件流线。 */
'use strict';
(function(){
 const data=window.ROCKET_GEOMETRY_CASES,query=new URLSearchParams(location.search),index=Number(query.get('case')||0),flow=query.get('flow')==='1';
 if(!Number.isInteger(index)||index<0||index>=data.length)throw Error('模型编号无效。');
 const select=document.getElementById('fixtureChoice');if(select)select.closest('label').hidden=true;
 const listeners=[];let value=null;
 const bridge={addEventListener:(name,fn)=>{if(name==='DataChanged')listeners.push(fn);}};
 Object.defineProperty(bridge,'Data',{get:()=>value,set:v=>{value=v;if(v&&['geometry-rendered','surface-flow-probe','geometry-error'].includes(v.kind))parent.postMessage({kind:'pressure-surface-view-event',caseIndex:index,payload:v},location.origin);}});
 window.setup(bridge);value={...data[index],showWire:query.get('wire')==='1'};
 if(!flow){delete value.surfaceFlowLines;value.coordinateNote='SI / m · 主轴X、鼻尖+X · 几何与网格，不表示物理解';}
 listeners.forEach(fn=>fn());
})();
