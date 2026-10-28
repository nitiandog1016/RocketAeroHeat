function benchmark=loadNasaBoattailBenchmark(folder)
%LOADNASABOATTAILBENCHMARK 核验NASA尾锥原图字节、读图标定及实验符号身份。
% 模块一：冻结来源清单与原PDF/图像/JSON，拒绝重绑或混入自造压力。
% 模块二：从原像素重算比例及5像素读图区间，实验总误差仍为未知。
% 模块三：图9横坐标与图1测孔表分别保存；不强行指定孔号，不生成外部CFD入口。
if nargin<1,folder=fullfile(fileparts(mfilename('fullpath')),'data','nasa_tnd2761');end
manifestPath=fullfile(folder,'source_manifest.json');
expected='9a523b5a4a882ff5488741a925fd7e71628699a46fc1d0c95b8627add0d47ddb';
assert(strcmp(fileSHA(manifestPath),expected),'NasaBoattail:Manifest','来源清单与已审查字节不同，禁止静默重绑。');
manifest=jsondecode(fileread(manifestPath));
for k=1:numel(manifest.files)
    f=manifest.files(k);name=char(f.path);
    assert(isempty(regexp(name,'(^[/\\])|(^[A-Za-z]:)|(^|[/\\])\.\.([/\\]|$)','once')), ...
        'NasaBoattail:Path','来源路径必须在本数据目录内。');
    p=fullfile(folder,name);info=dir(p);
    assert(~isempty(info)&&info.bytes==f.bytes&&strcmp(fileSHA(p),f.sha256), ...
        'NasaBoattail:Source','来源文件缺失或改字节：%s',name);
end
d=jsondecode(fileread(fullfile(folder,'pressure_digitization.json')));
g=jsondecode(fileread(fullfile(folder,'source_geometry.json')));
assert(strcmp(d.schema,'RocketAeroHeat.nasa-tnd2761.pressure-digitization.v1')&& ...
    d.alphaDeg==0&&d.noseSection==1&&d.pointCount==15&&numel(d.points)==15&& ...
    ~d.forcedOrificeMatch&&~d.completeNSRANSRunInputsAvailable&&~d.heatValidationComplete, ...
    'NasaBoattail:Scope','只能读取零攻角鼻型1的15个原图压力符号。');
assert(strcmp(d.sourcePDFSHA256,'8cc36817a2dcfb95a138e4436a0455bceaae02fbd7fe39c74c5a1bf0eef2c4d5')&& ...
    unknown(d.experimentalWallTemperatureK)&&unknown(d.fullExperimentalUncertainty), ...
    'NasaBoattail:Unknown','不能把缺失壁温或实验不确定性补成已知。');

%% 模块二：标定重算、未参与拟合网格检验与读图区间
c=d.calibration;A=[c.cornerPixels,ones(4,1)];matrix=A\c.cornerValues;
assert(max(abs(matrix-c.matrix),[],'all')<2e-13,'NasaBoattail:Calibration','像素标定矩阵不符。');
grid=[c.unusedGridPixels,ones(size(c.unusedGridPixels,1),1)]*matrix-c.unusedGridValues;
assert(max(abs(grid(:,1)))<.0032&&max(abs(grid(:,2)))<.004);
halfWidth=sum(abs(matrix(1:2,:)),1)*5;
assert(c.pixelReadAndRegistrationHalfWidth==5);
points=d.points;D=g.geometry.diameterM;Lbt=g.geometry.boattailAxialLengthM;
assert(abs(D-3.375*.0254)<1e-15&&abs(Lbt-1.18*.0254)<1e-15&& ...
    abs(g.geometry.assembledNoseBaseToBaseLengthM-14*.0254)<1e-15);
for k=1:numel(points)
    row=points(k);v=[row.pixelXY(:)' 1]*matrix;
    assert(row.noseSection==1&&row.alphaDeg==0&&ismember(row.boattailAngleDeg,[6 12 18])&& ...
        unknown(row.matchedOrificeID)&&unknown(row.absolutePressurePa)&&unknown(row.wallTemperatureK)&& ...
        ~row.digitizationIntervalIsStatisticalConfidence);
    assert(abs(v(1)-row.sourceFigureXprimeOverD)<2e-13&&abs(v(2)-row.pressureRatio)<2e-13);
    assert(abs(row.xDigitizationHalfWidth-halfWidth(1))<2e-13&&abs(row.ratioDigitizationHalfWidth-halfWidth(2))<2e-13);
    points(k).projectXFromFigureM=Lbt-D*v(1);
    points(k).projectXDigitizationHalfWidthM=D*halfWidth(1);
    points(k).coordinateIdentity='图9坐标转换；尚未与图1测孔逐点绑定';
    % jsondecode会把null变成空数组；显式恢复未知标量，重导出时保持null身份。
    for key={'matchedOrificeID','fullExperimentalUncertainty','absolutePressurePa','wallTemperatureK'}
        assert(unknown(row.(key{1})));points(k).(key{1})=NaN;
    end
end
for angle=[6 12 18]
    ids=find([points.boattailAngleDeg]==angle);assert(numel(ids)==5);
    assert(isequal([points(ids).sequenceAlongFigure],1:5)&& ...
        all(diff([points(ids).sourceFigureXprimeOverD])>0));
end
benchmark=struct('sourceManifestSHA256',expected,'sourceGeometry',g,'digitization',d, ...
    'points',points,'status','verified-source-calibration-and-experimental-symbol-transcription', ...
    'comparisonScope','无量纲压力关联对实验图符号的条件比较；非完整NS/RANS实验复现', ...
    'solverValidationComplete',false,'newProductionSolvedAreaM2',0);
end
function yes=unknown(value)
yes=isempty(value)||(isnumeric(value)&&isscalar(value)&&isnan(value));
end
function value=fileSHA(path)
fid=fopen(path,'rb');assert(fid>=0,'NasaBoattail:Missing','来源文件不存在。');guard=onCleanup(@()fclose(fid));b=fread(fid,Inf,'*uint8');clear guard;
md=java.security.MessageDigest.getInstance('SHA-256');md.update(b);value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
