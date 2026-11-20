/* 验收页读取真实图题、拾取结果、定量图片与下载表。该脚本只放入验证副本。 */
(() => {
  const node=id=>document.getElementById(id);
  const original=window.setup;
  async function inspect(event) {
    event.snapshotCaption=node('caseCondition').textContent;
    event.fieldTime=node('fieldTime').textContent;
    event.legend=[node('legendMin').textContent,node('legendMax').textContent];
    if(!event.snapshotCaption.includes('启动瞬态'))return;
    const scene=node('scene'),rect=scene.getBoundingClientRect();
    const pointer={clientX:rect.left+rect.width/2,clientY:rect.top+rect.height/2,button:0,pointerId:1};
    const capture=scene.setPointerCapture,click=HTMLAnchorElement.prototype.click;
    const figure=RocketExport.figurePlan;
    const createURL=URL.createObjectURL;
    let sourceBlob;
    try{
      scene.setPointerCapture=()=>{};
      scene.dispatchEvent(new PointerEvent('pointerdown',pointer));
      scene.dispatchEvent(new PointerEvent('pointerup',pointer));
      await Promise.resolve();
      event.probeTimes=Array.from(node('probe').querySelectorAll('.probeRow > small:first-of-type'),n=>n.textContent);
      event.probeFace=node('probeTitle').textContent;
      HTMLAnchorElement.prototype.click=()=>{};
      URL.createObjectURL=function(blob){sourceBlob=blob;return createURL.call(URL,blob);};
      RocketExport.figurePlan=(...args)=>{
        const plan=figure(...args);event.pngTitle=plan.texts[0].text;return plan;
      };
      node('save').click();
      event.pngName=node('saveDownload').download;
      event.csvName=node('saveSurfaceDownload').download;
      event.snapshotPNG=node('saveDownload').href;
      const csv=await sourceBlob.text();
      const rows=csv.trimEnd().split('\r\n').slice(1);
      event.surfaceTimeRows=rows.filter(row=>row.includes('"固定瞬态快照","0.00000780936086544745","'+event.fieldTime+'"')).length;
      node('saveSection').click();
      const link=node('saveSectionDownload');
      const section=decodeURIComponent(link.href.slice(link.href.indexOf(',')+1)).split('\r\n').slice(1);
      event.sectionTimeRows=section.filter(row=>row.includes('"'+event.fieldTime+'"')).length;
      event.sectionRows=section.length;
    }finally{
      scene.setPointerCapture=capture;HTMLAnchorElement.prototype.click=click;RocketExport.figurePlan=figure;
      URL.createObjectURL=createURL;
    }
  }
  window.setup=component=>{
    async function report(event){
      try{await inspect(event);}catch(error){event.observerError=String(error);}
      component.Data=event;
    }
    const offline=Boolean(window.ROCKET_CASE&&window.ROCKET_CASE.kind==='aerothermal-model');
    original({
      addEventListener:(...args)=>component.addEventListener(...args),
      get Data(){return offline?null:component.Data;},
      set Data(event){
        if(offline)return;
        if(event.kind==='rendered'){report(event);return;}
        component.Data=event;
      }
    });
    if(offline)requestAnimationFrame(()=>report({kind:'offline-proof'}));
  };
})();
