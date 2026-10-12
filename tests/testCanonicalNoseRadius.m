function testCanonicalNoseRadius()
% 模块用途：球头后接柱/锥时仍识别局部球头，并拒绝非球形鼻部。
addpath(fileparts(fileparts(mfilename('fullpath'))));
for kind={'sphere','hemisphereCylinder','bluntedCone'}
    [m,meta]=buildCanonicalBody(kind{1},struct('resolution',96));
    e=estimateNoseRadius(m);
    assert(e.available,'%s 球头未识别',kind{1});
    assert(abs(e.radius/meta.noseRadius-1)<1e-6);
    assert(e.fitSectionCount>=6 && e.outerDepthRadiusFraction>=0.08);
    % 非球形椭球的多尺度半径应变化，不能当球头。
    bad=m;bad.V(:,2)=1.3*bad.V(:,2);
    eBad=estimateNoseRadius(bad);assert(~eBad.available);
end
for kind={'cone','coneCylinder'}
    m=buildCanonicalBody(kind{1},struct('resolution',96));
    e=estimateNoseRadius(m);assert(~e.available,'尖锥被误识别为球头');
end
% 针对原缺陷：半球 12 圈以后为圆柱，原固定 24 圈必然跨越接缝。
[m,meta]=buildCanonicalBody('hemisphereCylinder',struct('resolution',48));
e=estimateNoseRadius(m);assert(e.available && abs(e.radius/meta.noseRadius-1)<1e-6);
fprintf('Canonical local spherical nose detection and non-spherical rejection passed.\n');
end
