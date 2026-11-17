/* 模块职责：原面返回研究态及热流来源展示；不修改共享渲染器或提升物理资格。 */
'use strict';
function applyKinsolCoupled(model,data){
  if(model.faces.length!==19752||data.originalFaces!==19752||data.identity.physicalQualification!==false)throw new Error('RESEARCH_IDENTITY');
  for(const field of data.fields){
    if(field.values.length!==19752||!field.values.every(Number.isFinite))throw new Error('INCOMPLETE_RESEARCH_FIELD');
    model.fields.push({...field,status:Array(19752).fill(1),methodIndex:Array(19752).fill(0),methodNames:[field.method],
      missingReason:Array(19752).fill(''),availability:'available',computedStateLabel:'未收敛返回研究态',
      scope:data.scope,stateNote:data.caseNote});
  }
  model.defaultField='kinsolReturnedQ';model.caseNote=data.caseNote;
  model.resultCondition={caption:'M5 · 30 km · α20° · 给定Tw300 K · KINSOL控制 · 未收敛'};
  model.researchIdentity={...model.researchIdentity,kinsolCoupled:data.identity};return model;
}
if(typeof module!=='undefined')module.exports={applyKinsolCoupled};
if(typeof window!=='undefined')applyKinsolCoupled(window.ROCKET_CASE,window.KINSOL_COUPLED);
