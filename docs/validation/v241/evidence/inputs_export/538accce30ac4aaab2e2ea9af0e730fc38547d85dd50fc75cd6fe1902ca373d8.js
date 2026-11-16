/* 模块职责：核对真正用于前端的研究身份、曲线共同范围和完整CSV；不声称浏览器验收。 */
'use strict';const fs=require('fs'),assert=require('assert'),path=require('path');
const api=require('../web/streamtube-research-v241.js'),d=JSON.parse(fs.readFileSync(process.argv[2],'utf8'));
api.validate(d);let values=0,lines=0;
for(let i=0;i<4;i++){const text=api.csv(d,i).trim().split('\n'),columns=text[0].slice(1,-1).split('","');assert.equal(text.length,d.tubes[i].rows.length+1);text.slice(1).forEach((line,j)=>{const row=line.slice(1,-1).split('","');assert.equal(row[1],'unqualified-research');assert.equal(row[row.length-1],'unknown');d.columns.forEach(k=>{const v=d.tubes[i].rows[j][k],s=row[columns.indexOf(k)];if(v===null)assert.equal(s,'');else{assert.equal(Number(s),v);values++;}});});lines+=text.length-1;}
for(const k of Object.keys(api.fields)){const r=api.range(d,k,0,true);d.tubes.forEach(t=>t.rows.forEach(s=>{if(s[k]!==null)assert(s[k]>=r[0]&&s[k]<=r[1]);}));}
const bad={...d,edgeAccepted:true};assert.throws(()=>api.validate(bad));assert.throws(()=>api.validate({...d,physicalQualification:true}));
const invalid=JSON.parse(JSON.stringify(d));invalid.tubes[0].rows[0].metricValid=0;assert.throws(()=>api.validate(invalid));
console.log(JSON.stringify({passed:true,tubes:4,CSVRows:lines,exactCSVNumbers:values,qualifiedIdentityRejections:2,invalidMetricRejected:true,browserInteractionTested:false}));
