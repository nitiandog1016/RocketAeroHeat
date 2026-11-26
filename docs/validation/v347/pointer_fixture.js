// 仅在验收页中发出按下与释放事件。坐标投影、前后遮挡和命中判定均由正式页面计算。
const productionSetup=window.setup;
window.setup=function(component){
  productionSetup(component);
  component.addEventListener('DataChanged',()=>{
    const message=component.Data;
    if(message?.kind!=='fixture-click')return;
    const canvas=document.getElementById('scene');
    const bounds=canvas.getBoundingClientRect();
    const event={clientX:bounds.left+message.x*bounds.width,clientY:bounds.top+message.y*bounds.height,
      shiftKey:false,button:0,pointerId:1};
    canvas.setPointerCapture=()=>{};
    canvas.onpointerdown(event);
    canvas.onpointerup(event);
    if(component.Data.kind==='fixture-click')component.Data={kind:'fixture-background'};
  });
};
