function packed=packRocketSteadyGraph(operator,options)
%PACKROCKETSTEADYGRAPH 将完整五量C=-Jgraph按控制体块CSR存储，保留所有原CV。
% 模块一：同线五带和全部侧向连接均纳入；同一控制体对的三角面贡献先求和。
% 模块二：仅作代数重排，五个守恒量始终一起移动，不修改物理网格/方程。
% 模块三：构造前核对持久存储及算法数组估算预算，不声称它等于进程峰值。
if nargin<2,options=struct;end
defaults=struct('ordering','symrcm','maximumPersistentBytes',768*1024^2,'maximumConstructionBytes',2048*1024^2);
assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(defaults))),'BlockGraph:Options','块图选项无效。');
for key=fieldnames(defaults)',if ~isfield(options,key{1}),options.(key{1})=defaults.(key{1});end,end
assert(ischar(options.ordering)&&ismember(options.ordering,{'natural','symrcm'}),'BlockGraph:Ordering','仅支持明确的自然或对称RCM代数排序。');
validateattributes(options.maximumPersistentBytes,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(options.maximumConstructionBytes,{'numeric'},{'real','finite','scalar','positive'});
assert(isfield(operator,'wlsNormalFactor')&&~isfield(operator,'steadyGraphPreconditionerReport'), ...
    'BlockGraph:Source','输入必须是带时间对角的已验证五带B，不能二次去正则。');
nf=operator.sourceFaces;nl=operator.layers;n=nf*nl;ns=numel(operator.sideMinus);assert(n<2^32&&nl>=3);
validateattributes(operator.pseudoTimeStepSeconds,{'numeric'},{'real','finite','scalar','positive'});
normal=(1:n-nf)';second=(1:n-2*nf)';self=(1:n)';
rawRows=[self;normal;normal+nf;second;second+2*nf;double(operator.sideMinus);double(operator.sidePlus)];
rawCols=[self;normal+nf;normal;second+2*nf;second;double(operator.sidePlus);double(operator.sideMinus)];
rawCount=numel(rawRows);pattern=spones(sparse(rawRows,rawCols,1,n,n));blockCount=nnz(pattern);
assert(issymmetric(pattern)&&all(diag(pattern)==1),'BlockGraph:Pattern','完整原面块连接应结构对称并有全部对角。');
persistentBytes=blockCount*(25*8+4)+(n+1)*4+n*4+n*25*8+2*n*8;
constructionBytes=2*blockCount*25*8+rawCount*64+blockCount*48+n*128+4096*25*8*4;
assert(persistentBytes<=options.maximumPersistentBytes&&constructionBytes<=options.maximumConstructionBytes, ...
    'BlockGraph:Budget','块图/因子或构造数组估算超预算，未分配五量块。');
if strcmp(options.ordering,'symrcm'),order=symrcm(pattern)';else,order=(1:n)';end
inverseOrder=zeros(n,1);inverseOrder(order)=1:n;
permuted=pattern(order,order);clear pattern rawRows rawCols;
[columns,owners]=find(permuted');clear permuted;
rowCounts=accumarray(owners,1,[n 1]);rowPointer=uint32([1;1+cumsum(rowCounts)]);
positionMap=sparse(columns,owners,(1:blockCount)',n,n);diagonalPosition=uint32(full(diag(positionMap)));
assert(all(double(diagonalPosition)>=double(rowPointer(1:n)))&&all(double(diagonalPosition)<double(rowPointer(2:n+1))));
values=zeros(5,5,blockCount);timer=tic;D=operator.diagonal;
for a=1:5,D(a,a,:)=D(a,a,:)-1/operator.pseudoTimeStepSeconds;end
add(self,self,D);clear D;
add(normal,normal+nf,operator.upper(:,:,normal));add(normal+nf,normal,operator.lower(:,:,normal+nf));
add(second,second+2*nf,operator.secondUpper(:,:,second));add(second+2*nf,second,operator.secondLower(:,:,second+2*nf));
add(double(operator.sideMinus),double(operator.sidePlus),operator.minusCoupling);
add(double(operator.sidePlus),double(operator.sideMinus),operator.plusCoupling);
assert(all(isfinite(values),'all'),'BlockGraph:Nonfinite','完整块图数值非有限。');
packed=struct('kind','complete-five-conservation-block-CSR-unshifted-steady-graph','cellCount',n,'sourceFaces',nf,'layers',nl, ...
    'rowPointer',rowPointer,'column',uint32(columns),'diagonalPosition',diagonalPosition,'values',values, ...
    'order',order,'inverseOrder',inverseOrder,'residualScale',operator.residualScale,'options',options, ...
    'packingReport',struct('rawBlockContributions',rawCount,'uniqueStoredBlocks',blockCount,'mergedDuplicateContributions',rawCount-blockCount, ...
        'minimumBlocksPerRow',min(rowCounts),'maximumBlocksPerRow',max(rowCounts),'meanBlocksPerRow',mean(rowCounts), ...
        'persistentFactorStorageUpperBoundBytes',persistentBytes,'algorithmArrayConstructionEstimateBytes',constructionBytes, ...
        'sourceOperatorAndLibraryWorkspaceExcluded',true,'processPeakMemoryClaimed',false,'packingSeconds',toc(timer), ...
        'allOriginalCellsAndFiveEquationsKept',true,'physicalEquationChanged',false,'newActualNSEvaluations',0,'newProductionSolvedAreaM2',0));
fprintf('BLOCK_GRAPH cells=%d raw=%d unique=%d persistentBound=%.0f constructionEstimate=%.0f\n',n,rawCount,blockCount,persistentBytes,constructionBytes);
    function add(physicalRows,physicalCols,pages)
        % 独立的局部变量名，避免嵌套函数共享批次数造成混合/尾批越界。
        for gFirst=1:4096:numel(physicalRows)
            gIDs=gFirst:min(gFirst+4095,numel(physicalRows));
            gIndex=sub2ind([n n],inverseOrder(physicalCols(gIDs)),inverseOrder(physicalRows(gIDs)));
            gPositions=full(positionMap(gIndex));assert(all(gPositions>0));
            [gUnique,~,gGroup]=unique(gPositions);
            if numel(gUnique)==numel(gPositions)
                values(:,:,gPositions)=values(:,:,gPositions)+pages(:,:,gIDs);
            else
                for gA=1:5,for gB=1:5
                    gSum=accumarray(gGroup,reshape(pages(gA,gB,gIDs),[],1),[numel(gUnique) 1]);
                    values(gA,gB,gUnique)=values(gA,gB,gUnique)+reshape(gSum,1,1,[]);
                end,end
            end
        end
    end
end
