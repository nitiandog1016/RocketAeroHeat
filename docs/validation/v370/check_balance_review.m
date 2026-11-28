function check_balance_review(out)
%CHECK_BALANCE_REVIEW 核对全域收支、六区热壁项及工况切换在正常侧栏的显示。
cfg=createCaseConfig('simpleLauncher');
app=rocketAeroHeatWorkbench(cfg);
guard=onCleanup(@()delete(app)); %#ok<NASGU>
app.Visible='off';
press(app.UserData.controls.load);press(app.UserData.controls.confirm);
press(app.UserData.controls.generate);press(app.UserData.workflowButtons(3));
press(app.UserData.controls.confirm);press(app.UserData.controls.openFlightRecords);
path=fullfile(cfg.projectRoot,'docs','validation','v370','evidence','balance.json');
b=jsondecode(fileread(path));
note=@()strjoin(app.UserData.controls.flightNote.Value,newline);
total=sprintf('整域边界与储能之和 %.2f kW',b.nativeBoundaryPlusStorage(2)/1000);
assert(contains(note(),total) && contains(note(),sprintf('%+.1f W',b.nativeExcessImbalanceW)));
% 切换区域时使用原生弱壁项，保留原 STL 面功率的独立读数。
sel=app.UserData.controls.flightRegion;
keys=fieldnames(b.regions);
for j=1:numel(keys)
    key=keys{j};sel.Value=key;sel.ValueChangedFcn(sel,[]);drawnow;
    assert(contains(note(),sprintf('求解器热壁边界项 %.3f W',b.regions.(key).nativeWeakWallW)));
    assert(contains(note(),sprintf('入壁功率 %.3f W',b.regions.(key).originalSurfaceW)) && ~contains(note(),total));
end
sel.Value='all';sel.ValueChangedFcn(sel,[]);drawnow;
mach=app.UserData.controls.flightRecord;
for m=[.85 1 1.2 2 .3]
    mach.Value=m;mach.ValueChangedFcn(mach,[]);drawnow;
    assert(contains(note(),total)==(m==.3));
end
s=app.UserData.workflowState;
assert(isempty(s.result) && ~s.steps(5).available);
proof=struct('sourceSHA256',cm011FileSHA256(path),'checkedRegions',{keys}, ...
    'checkedMach',[.3 .85 1 1.2 2],'finalResultAvailable',false);
writelines(string(jsonencode(proof)),fullfile(out,'balance.json'),'Encoding','UTF-8');
end

function press(button)
button.ButtonPushedFcn(button,[]);drawnow;
end
