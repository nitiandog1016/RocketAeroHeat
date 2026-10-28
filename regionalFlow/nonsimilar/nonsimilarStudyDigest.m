function hash=nonsimilarStudyDigest(study)
%NONSIMILARSTUDYDIGEST 绑定已重放的条件、站值和证据，防止映射时误用旧缓存。
% 这是可复核的数据一致性合同，不是对恶意修改的安全认证。
names={'kind','schemaVersion','freestream','noseRadius','gas','uniformWallTemperature', ...
    'inletAngleDeg','maximumAngularStepDeg','acceptedIntervalDeg','station','source', ...
    'physicalSolution','eligibleForSurfaceCoverage','experimentalValidation', ...
    'historicalPreRunSourceBinding','currentSourceReplayPassed','currentSourceHashes'};
bound=struct;
for k=1:numel(names)
    if ~isfield(study,names{k}),error('RocketAeroHeat:NonSimilarStudy:Binding','缺少站表绑定字段。');end
    bound.(names{k})=study.(names{k});
end
bytes=unicode2native(jsonencode(bound),'UTF-8');
digest=java.security.MessageDigest.getInstance('SHA-256');digest.update(bytes);
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2)',1,[]));
end
