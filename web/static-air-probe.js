/* 从已核对的节点帧取值。界面和 CSV 共享同一组带身份、单位和来源的读数。 */
(function () {
  'use strict';
  const sides = ['x_minus','x_plus','y_minus','y_plus','z_minus','z_plus'];
  const axes = ['x','y','z'];
  const grid = [0,.25,.5,.75,1];
  const missingMethod = '未记录';
  function check(ok, message) { if (!ok) throw Error(message); }
  function number(text) {
    check(typeof text === 'string' && text.trim() !== '' && Number.isFinite(Number(text)), '节点数值缺失或无效');
    return Number(text);
  }
  function walls(config) {
    const match = /^\s*MARKER_ISOTHERMAL\s*=\s*\(([^)]+)\)/m.exec(config.text);
    check(match, '缺少等温壁边界记录');
    const parts = match[1].split(',').map(x => x.trim()), result = new Map();
    for (let i = 0; i < parts.length; i += 2) {
      check(sides.includes(parts[i]) && !result.has(parts[i]), '壁面标记无法识别或重复');
      result.set(parts[i], number(parts[i+1]));
    }
    check(result.size === sides.length, '壁面记录不完整');
    return result;
  }
  function readFrame(file, worker) {
    // 归档格式为无引号的数值列和可选的方法名。原编号与文件行号分别处理。
    const [header, ...lines] = file.text.trim().split(/\r?\n/);
    const keys = header.split(',');
    const required = ['id', ...axes, 'T_K', 'p_Pa'];
    check(required.every(k => keys.includes(k)) && new Set(keys).size === keys.length, '缺少或重复节点列');
    const nodes = new Map(), positions = new Set();
    for (const line of lines) {
      const cells = line.split(',');
      check(cells.length === keys.length, '节点行长度有误');
      const cell = key => cells[keys.indexOf(key)];
      const id = number(cell('id')), xyz = axes.map(k => number(cell(k)));
      check(Number.isInteger(id) && id >= 0 && id < 125 && !nodes.has(id), '原节点编号无效或重复');
      check(xyz.every(x => grid.includes(x)) && !positions.has(xyz.join(',')), '节点坐标无效或重复');
      const method = cell('method') || '';
      check(!method || !worker || method === worker, '节点帧与工作对象记录的方法不一致');
      nodes.set(id, {id, xyz, temperature:number(cell('T_K')), pressure:number(cell('p_Pa')), method});
      positions.add(xyz.join(','));
    }
    check(nodes.size === 125, '节点帧不完整');
    return {nodes, source:file};
  }
  function openArchive(data) {
    const result = {};
    let reference;
    for (const key of ['default','stable']) {
      const raw = data[key], receipt = JSON.parse(raw.initialization.text);
      check(receipt.nodes === 125 && receipt.timeUpdates === 0 && receipt.nonlinearUpdates === 0,
            '来源并非静止空气初始化回读');
      const worker = receipt.workerMethod || '';
      const run = {wall:walls(raw.case), config:raw.case, receipt:raw.initialization, worker};
      for (const frame of ['initial','terminal']) {
        run[frame] = readFrame(raw[frame], worker);
        const nodes = run[frame].nodes;
        if (reference) {
          for (const [id, node] of nodes) {
            check(node.xyz.every((v,j) => v === reference.get(id).xyz[j]), '切换记录时相同编号的坐标发生变化');
          }
        } else reference = nodes;
      }
      result[key] = run;
    }
    return result;
  }
  function probe(archive, sample, frame, id) {
    const run = archive[sample], snapshot = run[frame], node = snapshot.nodes.get(id);
    check(node, '原节点编号不存在');
    const tags = sides.filter((_,i) => node.xyz[Math.floor(i/2)] === i%2);
    const boundary = tags.length > 0;
    const prescribed = new Set(tags.map(tag => run.wall.get(tag)));
    check(prescribed.size <= 1, '相交壁面的给定温度不唯一');
    const wallStatus = boundary ? '未知' : '不适用';
    const fields = [
      ['流体节点温度', node.temperature, 'K', '初始化回读', node.method || missingMethod, snapshot.source],
      ['流体节点压力', node.pressure, 'Pa', '初始化回读', '原生节点状态回读', snapshot.source],
      ['给定壁温', boundary ? run.wall.get(tags[0]) : null, 'K', boundary ? '给定' : '不适用',
       boundary ? 'MARKER_ISOTHERMAL' : '', boundary ? run.config : null],
      ['壁面热流', null, 'W/m²', wallStatus, '', null],
      ['求解固体壁温', null, 'K', wallStatus, '', null]
    ];
    return {sample, frame, node, tags, fields, snapshot:snapshot.source, receipt:run.receipt, worker:run.worker};
  }
  function csv(selection) {
    const p = selection;
    const rows = [['配置','源帧','原节点编号','X_m','Y_m','Z_m','边界标记','物理量','数值','单位','身份','方法',
                   '变量来源','变量来源_SHA256','节点帧','节点帧_SHA256','工作对象方法','工作对象记录','工作对象记录_SHA256']];
    for (const [label, value, unit, status, method, source] of p.fields) {
      rows.push([p.sample,p.frame,p.node.id,...p.node.xyz,p.tags.join('|'),label,value,unit,status,method,
                 source?.path,source?.sha256,p.snapshot.path,p.snapshot.sha256,p.worker || missingMethod,
                 p.receipt.path,p.receipt.sha256]);
    }
    const quote = v => '"' + String(v ?? '').replace(/"/g,'""') + '"';
    return '\uFEFF' + rows.map(row => row.map(quote).join(',')).join('\r\n') + '\r\n';
  }
  if (typeof module === 'object' && module.exports) module.exports = {openArchive,probe,csv};
  if (typeof document === 'undefined') return;

  // 同一个选择状态驱动 XY 截面、读数、文件身份及下载。切换档案保留原节点编号。
  const $ = id => document.getElementById(id);
  const el = (name, text) => { const x = document.createElement(name); x.textContent = text; return x; };
  let archive, selection, bridge;
  function deliver(kind) {
    const p = selection;
    if (bridge && p) bridge.Data = {kind, sample:p.sample, frame:p.frame, node:p.node.id,
      filename:'static-air-'+p.sample+'-'+p.frame+'-node-'+p.node.id+'.csv', csv:csv(p)};
  }
  function select(sample, frame, id) {
    const p = probe(archive,sample,frame,id);
    $('case').value = sample; $('frame').value = frame; $('node').value = id; $('z').value = p.node.xyz[2];
    $('probe').textContent = '原节点 ' + id;
    $('position').textContent = 'XYZ / m　' + p.node.xyz.join('，') + '　·　' + (p.tags.join(' / ') || '内部节点');
    $('values').replaceChildren();
    for (const [label,value,unit,status] of p.fields) {
      const row = document.createElement('tr');
      for (const v of [label,value ?? '—',unit,status]) row.append(el('td',v));
      $('values').append(row);
    }
    $('method').textContent = '温度反算方法　' + (p.node.method || missingMethod);
    $('sources').replaceChildren();
    for (const [role,file] of [['节点帧',p.snapshot],['工作对象 · '+(p.worker || missingMethod),p.receipt],
                               ['壁面配置',archive[sample].config]]) {
      $('sources').append(el('div',role+'　'+file.path+'　SHA256 '+file.sha256));
    }
    $('slice').replaceChildren();
    const nodes = Array.from(archive[sample][frame].nodes.values()).filter(n => n.xyz[2] === p.node.xyz[2]);
    for (const y of grid.slice().reverse()) {
      $('slice').append(el('span',y));
      for (const n of nodes.filter(n => n.xyz[1] === y).sort((a,b) => a.xyz[0]-b.xyz[0])) {
        const button = el('button', n.id);
        button.type='button'; button.setAttribute('aria-label','原节点 '+n.id+'，X '+n.xyz[0]+'，Y '+y);
        button.setAttribute('aria-pressed',String(n.id === id));
        button.onclick = () => select(sample,frame,n.id);
        $('slice').append(button);
      }
    }
    $('slice').append(el('span','Y / X'), ...grid.map(x => el('span',x)));
    selection = p;
    deliver('static-air-selection');
  }
  function changePlane() {
    const p = selection;
    const next = Array.from(archive[p.sample][p.frame].nodes.values()).find(n =>
      n.xyz[0] === p.node.xyz[0] && n.xyz[1] === p.node.xyz[1] && n.xyz[2] === Number($('z').value));
    select(p.sample,p.frame,next.id);
  }
  try {
    archive = openArchive(JSON.parse($('archive').textContent));
    grid.forEach(z => {const option = el('option',z); option.value=z; $('z').append(option);});
    select('default','initial',0);
    $('case').onchange = () => select($('case').value,selection.frame,selection.node.id);
    $('frame').onchange = () => select(selection.sample,$('frame').value,selection.node.id);
    $('node').onchange = () => {
      const id = $('node').valueAsNumber;
      if (!Number.isInteger(id) || id<0 || id>124) { $('node').value = selection.node.id; return; }
      select(selection.sample,selection.frame,id);
    };
    $('z').onchange = changePlane;
    $('download').onclick = () => {
      if (bridge) { deliver('static-air-export'); return; }
      const p = selection, a = document.createElement('a');
      const url = URL.createObjectURL(new Blob([csv(p)], {type:'text/csv;charset=utf-8'}));
      a.download = 'static-air-'+p.sample+'-'+p.frame+'-node-'+p.node.id+'.csv'; a.href=url;
      document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url),1000);
    };
  } catch (err) {
    $('error').textContent = err.message; $('error').hidden = false;
    for (const id of ['case','frame','node','z','download']) $(id).disabled = true;
  }
  // MATLAB 可同步选择节点和源帧；单独打开时使用相同的选择、渲染和 CSV 路径。
  window.setup = component => {
    bridge=component;
    component.addEventListener('DataChanged',() => {
      const m = component.Data;
      if (m?.kind !== 'static-air-select' || !selection) return;
      const sample=m.sample ?? selection.sample, frame=m.frame ?? selection.frame, id=m.node ?? selection.node.id;
      if (['default','stable'].includes(sample) && ['initial','terminal'].includes(frame)
          && Number.isInteger(id) && id>=0 && id<125) select(sample,frame,id);
    });
    deliver('static-air-selection');
  };
})();
