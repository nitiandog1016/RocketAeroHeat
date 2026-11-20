/* 仅随验证窗口加载。读取实际色标与 PNG 导出结果，文件保存在验证目录。 */
(() => {
  const setupViewer = window.setup;
  window.setup = component => {
    setupViewer({
      addEventListener: (...args) => component.addEventListener(...args),
      get Data() { return component.Data; },
      set Data(event) {
        if (event.kind === 'rendered') {
          const get = id => document.getElementById(id);
          event.legend = {low:get('legendMin').textContent,high:get('legendMax').textContent,
            note:get('legendNote').textContent,swatch:document.querySelector('.gradient').style.background};
          const caption=get('caseCondition').textContent;
          if (caption.includes('给定壁温 293.15 K') || caption.includes('给定壁温 299.15 K')) {
            const plan = RocketExport.figurePlan;
            const click = HTMLAnchorElement.prototype.click;
            try {
              // 执行原 PNG 和来源 CSV 回调，只截住最后的系统下载动作。
              HTMLAnchorElement.prototype.click = function () {};
              RocketExport.figurePlan = function (model,field,options,measure) {
                const result=plan(model,field,options,measure);
                event.exportRange=result.range;event.exportColors=result.colors;
                return result;
              };
              get('save').click();
              event.png=get('saveDownload').href;
            } finally {
              RocketExport.figurePlan=plan;HTMLAnchorElement.prototype.click=click;
            }
          }
        }
        component.Data=event;
      }
    });
  };
})();
