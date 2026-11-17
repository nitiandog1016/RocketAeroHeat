/* 模块职责：完整域零时间接口诊断接原面，不把研究压力或残差当合格气动热。 */
'use strict';
function applyFullResidual(model,data){
  if(model.faces.length!==19752||data.originalFaces!==19752)throw new Error('ORIGINAL_FACE_COUNT');
  for(const field of data.fields){
    if(field.values.length!==19752||!field.values.every(Number.isFinite))throw new Error('INCOMPLETE_DIAGNOSTIC');
    model.fields.push({...field,status:Array(19752).fill(1),methodIndex:Array(19752).fill(0),methodNames:[field.method],
      missingReason:Array(19752).fill(''),availability:'available',computedStateLabel:'未收敛研究态／算子诊断',
      scope:data.scope,stateNote:'零时间评估，无非线性推进；数据覆盖不等于物理合格覆盖。'});
  }
  model.defaultField='nativePhysicalEnergyResidual';model.caseNote=data.caseNote;
  model.resultCondition={caption:'M5 · 30 km · α20° · 给定Tw300 K · 原完整域 · 零时间接口检查'};
  model.researchIdentity={...model.researchIdentity,fullResidual:data.identity};
  return model;
}
if(typeof module!=='undefined')module.exports={applyFullResidual};
if(typeof window!=='undefined')applyFullResidual(window.ROCKET_CASE,window.FULL_RESIDUAL);
