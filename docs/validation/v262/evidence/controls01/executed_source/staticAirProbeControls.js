/* 仅在测试副本中接收真实 DOM 控件的状态及下载字节，不改变交付页面。 */
window.setup = async component => {
  const get = id => document.getElementById(id);
  const change = (id,value) => {
    const input=get(id);input.value=value;input.dispatchEvent(new Event('change',{bubbles:true}));
  };
  try {
    change('z','0.5');
    const target=get('slice').querySelector('button[aria-label^="原节点 62，"]');
    if (!target) throw Error('缺少内部节点按钮');
    target.click();
    change('case','stable');change('frame','terminal');
    const table=Array.from(get('values').rows,r=>Array.from(r.cells,c=>c.textContent));
    const name=get('probe').textContent, method=get('method').textContent, position=get('position').textContent;
    const selected=get('slice').querySelector('button[aria-pressed="true"]').textContent;
    const sources=get('sources').textContent;
    const download=new Promise((resolve,reject)=>{
      document.addEventListener('click', async event=>{
        if (!(event.target instanceof HTMLAnchorElement)) return;
        event.preventDefault();
        try {
          const response=await fetch(event.target.href);
          const bytes=Array.from(new Uint8Array(await response.arrayBuffer()));
          resolve({filename:event.target.download,bytes});
        } catch (err) { reject(err); }
      },true);
    });
    get('download').click();
    const file=await download;
    component.Data={ok:true,cells:JSON.stringify(table.flat()),name,method,position,selected,sources,
                    filename:file.filename,bytes:file.bytes};
  } catch (err) {component.Data={ok:false,error:err.message};}
};
