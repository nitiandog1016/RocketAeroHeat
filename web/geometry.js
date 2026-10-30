/* RocketAeroHeat 离线几何预览；可选逐面研究场量，未知值保持灰色和可拾取。 */
'use strict';
(function(){
  const $=id=>document.getElementById(id),canvas=$('scene'),finite=x=>typeof x==='number'&&Number.isFinite(x);
  const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2],fmt=x=>Number(x.toPrecision(5)).toString();
  const fieldColor=t=>{const stops=[[.94,.97,.99],[.67,.82,.90],[.28,.60,.75],[.08,.35,.56],[.025,.15,.31]],x=Math.max(0,Math.min(1,t))*4,i=Math.min(3,Math.floor(x)),f=x-i;return stops[i].map((v,j)=>v+(stops[i+1][j]-v)*f);};
  let gl,program,positions,colorBuffer,edges,flowPositions,flowColors,sectionPositions,bridge,model=null,vertices=[],lo,hi,initializationError='';
  let modelCenter=[0,0,0],modelExtent=1,flowSegments=0,flowPoints=[],sectionSegmentCount=0;
  let az=-Math.PI/2,el=.35,zoom=1,pan=[0,0],baseScale=1,right,up,direction,sx,sy,frame=0,pendingReport=false,drag;
  // 模块一：独立 WebGL 管线，仅使用几何法向和固定中性色，不接收物理字段。
  try{
    gl=canvas.getContext('webgl',{antialias:true,preserveDrawingBuffer:true});
    if(!gl)throw Error('当前内嵌网页环境不支持 WebGL。');
    const shader=(type,source)=>{const s=gl.createShader(type);gl.shaderSource(s,source);gl.compileShader(s);if(!gl.getShaderParameter(s,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(s));return s;};
    program=gl.createProgram();
    gl.attachShader(program,shader(gl.VERTEX_SHADER,'attribute vec3 p,c;uniform vec3 r,u,d;uniform vec2 scale,pan;varying vec3 color;void main(){gl_Position=vec4(dot(p,r)*scale.x+pan.x,dot(p,u)*scale.y+pan.y,-dot(p,d)*0.45,1.0);color=c;}'));
    gl.attachShader(program,shader(gl.FRAGMENT_SHADER,'precision mediump float;varying vec3 color;void main(){gl_FragColor=vec4(color,1.0);}'));
    gl.linkProgram(program);if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(program));
    positions=gl.createBuffer();colorBuffer=gl.createBuffer();edges=gl.createBuffer();flowPositions=gl.createBuffer();flowColors=gl.createBuffer();sectionPositions=gl.createBuffer();gl.enable(gl.DEPTH_TEST);gl.clearColor(0,0,0,0);
  }catch(e){initializationError=e.message;problem(e.message);}
  function problem(text){$('error').hidden=false;$('error').textContent=text;if(bridge)bridge.Data={kind:'geometry-error',revision:model?model.revision:0,message:text,nonce:Date.now()};}
  const buffer=(b,v)=>{gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.bufferData(gl.ARRAY_BUFFER,new Float32Array(v),gl.STATIC_DRAW);};
  function attribute(name,b){const i=gl.getAttribLocation(program,name);gl.bindBuffer(gl.ARRAY_BUFFER,b);gl.enableVertexAttribArray(i);gl.vertexAttribPointer(i,3,gl.FLOAT,false,0,0);}
  // 模块用途：可选的条件流线叠加；几何仍保持中性色，不把输入p/T绘成已解物理场。
  function loadFlow(overlay){
    flowSegments=0;flowPoints=[];const ps=[],cs=[];
    if(overlay){
      if(overlay.kind!=='pressure-surface-flow-research'||overlay.productionQualified!==false||!Array.isArray(overlay.curves))throw Error('表面流线研究身份无效。');
      overlay.curves.forEach((curve,curveIndex)=>{
        const pts=curve.points,ns=curve.normals;
        if(!Array.isArray(pts)||!Array.isArray(ns)||pts.length!==ns.length||pts.some(p=>!Array.isArray(p)||p.length!==3||p.some(x=>!finite(x)))||ns.some(p=>!Array.isArray(p)||p.length!==3||p.some(x=>!finite(x))||Math.hypot(...p)===0))throw Error('流线坐标或法向无效。');
        const normalized=pts.map(p=>p.map((x,j)=>(x-modelCenter[j])/modelExtent));
        normalized.forEach((p,stationIndex)=>flowPoints.push({p,normal:ns[stationIndex],curveIndex,stationIndex}));
        for(let k=1;k<normalized.length;k++){ps.push(...normalized[k-1],...normalized[k]);cs.push(.025,.25,.4,.025,.25,.4);flowSegments++;}
      });
    }
    if(gl){buffer(flowPositions,ps);buffer(flowColors,cs);}
    const note=document.querySelector('footer span:last-child');
    if(note)note.textContent=overlay?(overlay.thermodynamicStateIsConditionallyRestored===true?'条件流线数值 · 输入p / 同物性反算Te · Tr/Tw/热流未求':'条件流线数值 · 压力/热状态为输入 · 热流未求'):'几何预览 · 尚未计算物理场量';
  }
  function camera(){direction=[Math.cos(el)*Math.cos(az),Math.cos(el)*Math.sin(az),Math.sin(el)];right=[-Math.sin(az),Math.cos(az),0];up=[-Math.sin(el)*Math.cos(az),-Math.sin(el)*Math.sin(az),Math.cos(el)];sx=baseScale*zoom;sy=sx*canvas.clientWidth/Math.max(1,canvas.clientHeight);}
  function fit(){camera();if(!vertices.length)return;let w=0,h=0;vertices.forEach(v=>{w=Math.max(w,Math.abs(dot(v,right)));h=Math.max(h,Math.abs(dot(v,up)));});const aspect=canvas.clientWidth/Math.max(1,canvas.clientHeight);baseScale=Math.min(.82/Math.max(w,1e-8),.76/Math.max(h*aspect,1e-8));zoom=1;pan=[0,0];}
  // 模块二：真实三角面展开；统一缩放 XYZ，保持长度/直径的实际关系。
  function load(data){
    try{
      if(!data||data.kind!=='aerothermal-geometry'||data.schemaVersion!==1||!Number.isInteger(data.revision))throw Error('几何负载类型或版本无效。');
      let V=data.vertices,F=data.faces;
      if(Array.isArray(V)&&typeof V[0]==='number')V=[V];if(Array.isArray(F)&&typeof F[0]==='number')F=[F];
      if(!Array.isArray(V)||!Array.isArray(F)||Boolean(V.length)!==Boolean(F.length)||V.some(v=>!Array.isArray(v)||v.length!==3||v.some(x=>!finite(x)))||F.some(f=>!Array.isArray(f)||f.length!==3||f.some(i=>!Number.isInteger(i)||i<0||i>=V.length)))throw Error('顶点坐标或三角形编号无效。');
      const field=data.faceField;
      if(field&&(!Array.isArray(field.values)||field.values.length!==F.length||field.values.some(v=>v!==null&&!finite(v))||!Array.isArray(field.range)||field.range.length!==2||!field.range.every(finite)||field.range[1]<=field.range[0]||typeof field.name!=='string'||typeof field.unit!=='string'||field.productionQualified!==false))throw Error('逐面研究场量、单位、色标或资格无效。');
      model={vertices:V,faces:F,revision:data.revision,faceField:field||null,hasSectionOverlay:Array.isArray(data.sectionSegments)};
      if(initializationError)throw Error(initializationError);
      vertices=[];lo=[Infinity,Infinity,Infinity];hi=[-Infinity,-Infinity,-Infinity];
      V.forEach(v=>v.forEach((x,j)=>{lo[j]=Math.min(lo[j],x);hi[j]=Math.max(hi[j],x);}));
      if(V.length){const extent=Math.max(...hi.map((x,j)=>x-lo[j]));if(!(extent>0))throw Error('模型尺寸为零。');const center=lo.map((x,j)=>(x+hi[j])/2);modelCenter=center;modelExtent=extent;vertices=V.map(v=>v.map((x,j)=>(x-center[j])/extent));}
      loadFlow(data.surfaceFlowLines);
      // 模块用途：仅叠加真实空间截面交线，保持原相机和面元拾取，不把交线当流线。
      const sections=data.sectionSegments||[];
      if(!Array.isArray(sections)||sections.length>100000||sections.some(s=>!Array.isArray(s)||s.length!==2||s.some(p=>!Array.isArray(p)||p.length!==3||p.some(x=>!finite(x)))))throw Error('截面交线坐标无效。');
      sectionSegmentCount=sections.length;
      if(gl)buffer(sectionPositions,sections.flatMap(s=>s.flatMap(p=>p.map((x,j)=>(x-modelCenter[j])/modelExtent))));
      if(gl){const ts=[],cs=[],es=[];F.forEach((f,faceIndex)=>{const a=vertices[f[0]],b=vertices[f[1]],c=vertices[f[2]],ab=b.map((x,j)=>x-a[j]),ac=c.map((x,j)=>x-a[j]);let n=[ab[1]*ac[2]-ab[2]*ac[1],ab[2]*ac[0]-ab[0]*ac[2],ab[0]*ac[1]-ab[1]*ac[0]];const norm=Math.hypot(...n);n=norm>0?n.map(x=>x/norm):[0,0,0];const light=.65+.22*Math.abs(dot(n,[.3,.4,.8660254])),color=field?(finite(field.values[faceIndex])?fieldColor((field.values[faceIndex]-field.range[0])/(field.range[1]-field.range[0])):[.66,.69,.72]):[.52,.60,.67].map(x=>x*light);f.forEach(i=>{ts.push(...vertices[i]);cs.push(...color);});[0,1,1,2,2,0].forEach(j=>es.push(...vertices[f[j]]));});buffer(positions,ts);buffer(colorBuffer,cs);buffer(edges,es);}
      const legend=$('fieldLegend');if(legend){legend.hidden=!field;if(field){legend.replaceChildren();const label=document.createElement('span'),bar=document.createElement('div');label.textContent=field.name+' / '+field.unit+' · 固定范围 '+fmt(field.range[0])+'—'+fmt(field.range[1])+' · 超范围使用端点色；灰色为未知';bar.className='fieldColorBar';legend.append(label,bar);}}
      $('caption').textContent=data.caption||'几何预览';$('coordinates').textContent=data.coordinateNote||'SI · XYZ / m · 主轴 X · 鼻尖 +X · 真实几何比例';$('wire').checked=Boolean(data.showWire);
      $('summary').textContent=F.length?V.length.toLocaleString()+' 顶点 · '+F.length.toLocaleString()+' 三角面':'等待导入模型';$('empty').hidden=Boolean(F.length);
      $('bounds').textContent=V.length?'几何尺寸 / m：X '+fmt(hi[0]-lo[0])+' · Y '+fmt(hi[1]-lo[1])+' · Z '+fmt(hi[2]-lo[2]):'几何尺寸 / m：—';
      if(data.resetCamera){const view=finite(data.initialCameraAzimuthRad)?data:data.surfaceFlowLines;az=view&&finite(view.initialCameraAzimuthRad)?view.initialCameraAzimuthRad:-Math.PI/2;el=view&&finite(view.initialCameraElevationRad)?Math.max(-1.55,Math.min(1.55,view.initialCameraElevationRad)):.35;fit();}if(gl)$('error').hidden=true;pendingReport=true;render();
    }catch(e){problem(e.message);}
  }
  // 模块三：正交投影与像素回执；透明清屏便于核验实际非背景像素。
  function render(){if(!frame)frame=requestAnimationFrame(draw);}
  function draw(){frame=0;if(!gl||!model)return;const dpr=Math.min(window.devicePixelRatio||1,2),w=Math.max(1,Math.round(canvas.clientWidth*dpr)),h=Math.max(1,Math.round(canvas.clientHeight*dpr));if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}gl.viewport(0,0,w,h);gl.clear(gl.COLOR_BUFFER_BIT|gl.DEPTH_BUFFER_BIT);gl.useProgram(program);camera();
    [['r',right],['u',up],['d',direction]].forEach(([k,v])=>gl.uniform3fv(gl.getUniformLocation(program,k),v));gl.uniform2f(gl.getUniformLocation(program,'scale'),sx,sy);gl.uniform2fv(gl.getUniformLocation(program,'pan'),pan);
    attribute('p',positions);attribute('c',colorBuffer);gl.enable(gl.POLYGON_OFFSET_FILL);gl.polygonOffset(1,1);gl.drawArrays(gl.TRIANGLES,0,model.faces.length*3);gl.disable(gl.POLYGON_OFFSET_FILL);
    if($('wire').checked){attribute('p',edges);const i=gl.getAttribLocation(program,'c');gl.disableVertexAttribArray(i);gl.vertexAttrib3fv(i,[.17,.23,.29]);gl.drawArrays(gl.LINES,0,model.faces.length*6);}
    if(flowSegments){attribute('p',flowPositions);attribute('c',flowColors);gl.drawArrays(gl.LINES,0,flowSegments*2);}
    if(sectionSegmentCount){attribute('p',sectionPositions);const i=gl.getAttribLocation(program,'c');gl.disableVertexAttribArray(i);gl.vertexAttrib3fv(i,[.74,.29,.10]);gl.drawArrays(gl.LINES,0,sectionSegmentCount*2);}
    $('axes').replaceChildren();['X','Y','Z'].forEach((t,i)=>{const s=document.createElement('span');s.textContent=t;s.style.left=(32+right[i]*27)+'px';s.style.top=(32-up[i]*27)+'px';s.style.color=['#9c5046','#427357','#466d92'][i];$('axes').append(s);});
    if(pendingReport&&bridge){pendingReport=false;const pixels=new Uint8Array(w*h*4);gl.readPixels(0,0,w,h,gl.RGBA,gl.UNSIGNED_BYTE,pixels);let visible=0;for(let i=3;i<pixels.length;i+=64)if(pixels[i])visible++;bridge.Data={kind:'geometry-rendered',revision:model.revision,faces:model.faces.length,vertices:model.vertices.length,visiblePixelSamples:visible,webglError:gl.getError(),width:w,height:h,firstVertex:model.vertices[0]||[],bounds:model.vertices.length?[lo,hi]:[],...(flowSegments?{surfaceFlowSegments:flowSegments,surfaceFlowPointCount:flowPoints.length}:{}),...(model.hasSectionOverlay?{sectionSegments:sectionSegmentCount}:{}),nonce:Date.now()};}
  }
  // 模块四：旋转/平移/缩放及标准视角；无持续动画或外部资源请求。
  canvas.oncontextmenu=e=>e.preventDefault();canvas.onpointerdown=e=>{drag={x:e.clientX,y:e.clientY,startX:e.clientX,startY:e.clientY,pan:e.shiftKey||e.button===2||$('panMode').checked};canvas.setPointerCapture(e.pointerId);};
  canvas.onpointermove=e=>{if(!drag)return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;drag.x=e.clientX;drag.y=e.clientY;if(drag.pan){pan[0]+=2*dx/Math.max(1,canvas.clientWidth);pan[1]-=2*dy/Math.max(1,canvas.clientHeight);}else{az-=dx*.008;el=Math.max(-1.55,Math.min(1.55,el+dy*.008));}render();};
  canvas.onpointerup=e=>{
    const pick=drag&&!drag.pan&&Math.hypot(e.clientX-drag.startX,e.clientY-drag.startY)<4;drag=null;
    if(!pick||!bridge)return;camera();const rect=canvas.getBoundingClientRect();let chosen=null,best=12*12,depth=-Infinity;
    // 模块用途：按真实三角面投影重心及前后深度拾取，不要求场量有限，不让未知区消失。
    if(model&&model.faceField&&!flowPoints.length){let hit=null;model.faces.forEach((f,faceIndex)=>{const p=f.map(i=>vertices[i]),xy=p.map(v=>[rect.left+(1+dot(v,right)*sx+pan[0])*rect.width/2,rect.top+(1-dot(v,up)*sy-pan[1])*rect.height/2]);const a=xy[0],b=xy[1],c=xy[2],den=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]);if(Math.abs(den)<1e-12)return;const u=((b[1]-c[1])*(e.clientX-c[0])+(c[0]-b[0])*(e.clientY-c[1]))/den,v=((c[1]-a[1])*(e.clientX-c[0])+(a[0]-c[0])*(e.clientY-c[1]))/den,w=1-u-v;if(Math.min(u,v,w)<-1e-9)return;const bary=[u,v,w],z=bary.reduce((s,t,i)=>s+t*dot(p[i],direction),0);if(z>depth){depth=z;hit={faceIndex,barycentric:bary,position:[0,1,2].map(j=>bary.reduce((s,t,i)=>s+t*model.vertices[f[i]][j],0))};}});if(hit)bridge.Data={kind:'geometry-face-probe',revision:model.revision,...hit,nonce:Date.now()};return;}
    if(!flowPoints.length)return;
    flowPoints.forEach(point=>{if(dot(point.normal,direction)<=0)return;const x=rect.left+(1+dot(point.p,right)*sx+pan[0])*rect.width/2,y=rect.top+(1-dot(point.p,up)*sy-pan[1])*rect.height/2;const distance=(x-e.clientX)**2+(y-e.clientY)**2,z=dot(point.p,direction);if(distance<best||(distance===best&&z>depth)){chosen=point;best=distance;depth=z;}});
    if(chosen)bridge.Data={kind:'surface-flow-probe',revision:model.revision,curveIndex:chosen.curveIndex,stationIndex:chosen.stationIndex,nonce:Date.now()};
  };canvas.onpointercancel=()=>{drag=null;};canvas.addEventListener('wheel',e=>{e.preventDefault();zoom=Math.max(.05,Math.min(30,zoom*Math.exp(-e.deltaY*.001)));render();},{passive:false});
  $('wire').onchange=render;$('reset').onclick=()=>{fit();render();};document.querySelectorAll('[data-view]').forEach(b=>b.onclick=()=>{const v=b.dataset.view;az=v==='nose'?0:v==='iso'?-1.0:-Math.PI/2;el=v==='nose'||v==='side'?0:.5;fit();render();});
  const observer=new ResizeObserver(()=>{fit();render();});observer.observe(canvas);
  // 模块五：MATLAB uihtml 双向白名单桥接；收到确认事件不重复加载几何。
  window.setup=function(component){bridge=component;component.addEventListener('DataChanged',()=>{if(component.Data&&component.Data.kind==='aerothermal-geometry')load(component.Data);});component.Data={kind:'ready',nonce:Date.now()};};
})();
