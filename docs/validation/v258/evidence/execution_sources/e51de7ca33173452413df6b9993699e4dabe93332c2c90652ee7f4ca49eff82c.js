/* 模块职责：在原面显示固定方向诊断；差分未获一致性，不把它发布为精确Jacobian或新热流。 */
'use strict';
function applyDirectionConsistency(model,data){
  if(model.faces.length!==19752||data.originalFaces!==19752||data.identity.physicalQualification!==false)throw new Error('DIAGNOSTIC_IDENTITY');
  for(const field of data.fields){
    if(field.values.length!==19752||!field.values.every(Number.isFinite))throw new Error('INCOMPLETE_DIRECTION_FIELD');
    model.fields.push({...field,status:Array(19752).fill(1),methodIndex:Array(19752).fill(0),methodNames:[field.method],
      missingReason:Array(19752).fill(''),availability:'available',computedStateLabel:'方向诊断／差分一致性未通过',
      scope:data.scope,stateNote:data.caseNote});
  }
  model.defaultField='directionEnergyDifference';model.caseNote=data.caseNote;
  model.resultCondition={caption:'M5 · 30 km · α20° · 固定未接受方向诊断 · 无新流场更新'};
  model.researchIdentity={...model.researchIdentity,directionConsistency:data.identity};return model;
}
if(typeof module!=='undefined')module.exports={applyDirectionConsistency};
if(typeof window!=='undefined')applyDirectionConsistency(window.ROCKET_CASE,window.DIRECTION_CONSISTENCY);
