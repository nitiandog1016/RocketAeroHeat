function data=readSchaeferHeatTransferData(file)
%READSCHAEFERHEATTRANSFERDATA 读取原图11的成组热流/壁温草稿，保留像素、缺测和不确定度。
% 模块一：核对原图SHA、线性温度轴和对数St轴；同一圆/方系列分别配对。
% 模块二：像素横坐标独立核对图1名义测站；不移动原点或放宽容差制造通过。
% 模块三：输出St及测量温度，不自动换成q，不延拓未测鼻端/底面，不授予实验资格。
if nargin<1,file=fullfile(fileparts(mfilename('fullpath')),'tests','data','schaefer_ferguson_1961','heat_fig11.json');end
raw=jsondecode(fileread(file));root=fileparts(mfilename('fullpath'));
assert(raw.schemaVersion==1 && raw.figure==11 && ~raw.physicsQualified);
source=fullfile(root,raw.archivedSourceImage);assert(strcmp(sha(source),raw.sourceImageSHA256),'实验原图身份变化。');
c=raw.calibration;fit=polyfit(log10(c.heatStSqrtRe),c.heatY,1);
assert(fit(1)<0 && max(abs(polyval(fit,log10(c.heatStSqrtRe))-c.heatY))<1);
n=numel(raw.points);St=nan(n,1);Tw=St;Taw=St;X=St;plotX=St;delta=nan(n,3);unc=nan(n,2);
missing=false(n,1);broken=missing;ambiguous=missing;match=missing;ids=cell(n,1);series=ids;
T0=raw.condition.totalTemperature_R*5/9;D=raw.condition.diameter_in*.0254;L=17.57*.0254;
for k=1:n
    p=raw.points(k);ids{k}=p.id;series{k}=p.series;X(k)=L-p.station_in*.0254;
    missing(k)=isempty(p.heatPixel)||isempty(p.wallPixel);broken(k)=strcmp(p.symbolStatus,'broken-symbol');
    ambiguous(k)=strcmp(p.symbolStatus,'overlapping-symbol-or-line');
    Taw(k)=ratio(p.adiabaticPixel(2),c)*T0;
    expectedUpper=c.upperX(1)+diff(c.upperX)*(p.station_in/raw.condition.diameter_in)/10;
    expectedLower=c.lowerX(1)+diff(c.lowerX)*(p.station_in/raw.condition.diameter_in)/10;
    delta(k,2)=abs(expectedUpper-p.adiabaticPixel(1));
    if missing(k),continue;end
    Tw(k)=ratio(p.wallPixel(2),c)*T0;
    ordinate=10^((p.heatPixel(2)-fit(2))/fit(1));St(k)=ordinate/sqrt(raw.condition.diameterReynolds);
    % 对数坐标的像素误差转成非对称区间，只代表读图误差，不替代20%测量误差。
    u=raw.manualCenterUncertainty_pixels;
    unc(k,:)=sort(10.^(([p.heatPixel(2)-u p.heatPixel(2)+u]-fit(2))/fit(1))/sqrt(raw.condition.diameterReynolds));
    plotX(k)=L-(p.heatPixel(1)-c.lowerX(1))/diff(c.lowerX)*10*D;
    delta(k,1)=abs(expectedUpper-p.wallPixel(1));delta(k,3)=abs(expectedLower-p.heatPixel(1));
    match(k)=all(delta(k,:)<=raw.stationTolerance_pixels);
end
assert(all(Taw>0) && all(Tw(~missing)>0) && all(Taw(~missing)>Tw(~missing)));
tableOut=table(ids,series,X,plotX,Tw,Taw,St,unc(:,1),unc(:,2),missing,broken,ambiguous,match, ...
    'VariableNames',{'id','series','nominalX_m','heatPlotX_m','Tw_K','Taw_K','St', ...
    'St_digitizationLower','St_digitizationUpper','missing','broken','overlapNeedsReview','nominalStationMatched'});
data=struct('raw',raw,'points',tableOut,'stationDifferencePixels',delta,'logAxisFit',fit, ...
    'temperatureRatioDigitizationUncertainty',raw.manualCenterUncertainty_pixels/abs(diff(c.temperatureY)), ...
    'inputStatus','配对读图草稿；测站/重叠需复审，缺测保留，不能直接作全表面边界', ...
    'heatFluxConversionPerformed',false,'wallTemperatureExtended',false,'experimentalComparisonPerformed',false, ...
    'physicsQualified',false,'newQualifiedArea',0);
end
function t=ratio(y,c)
t=c.temperatureRatio(1)+(y-c.temperatureY(1))*diff(c.temperatureRatio)/diff(c.temperatureY);
end
function h=sha(file)
fid=fopen(file,'rb');assert(fid>=0);guard=onCleanup(@()fclose(fid));md=java.security.MessageDigest.getInstance('SHA-256');
while ~feof(fid),md.update(fread(fid,1024*1024,'*uint8'));end
h=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
